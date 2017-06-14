# Script to add componenet to the image and to do completeness test
# Copyleft (c) Ridlo W. Wibowo - 2017

import os
import sys
import glob
import numpy as np
#from astropy import units as u 
#from astropy import coordinates

class completenessTest:
    def __init__(self):
        pass

    
    def generate_sample(self, flux=1.02e-04, freq='241.23GHz', center=[65.8158364, -1.3425182], 
        PB=24.0, beam=0.2, imagesize=48.0, spacing=10, outputfile="sample.dat"):
        """
        Generate grid sample, with spacing parameter as input
        spacing is multiple of beamsize
        """

        PB_in_deg = PB/3600.0

        distance_between_source = spacing*beam

        grid = imagesize/(distance_between_source)
        
        number = (grid)**2
        
        print grid, number
        
        grid = int(grid)
        print "Grid size : ", grid, "x", grid, "\nTotal number of sample: ", grid*grid

        start  = [center[0] - imagesize/(2*3600) + 0.5*distance_between_source/3600.0, center[1] - imagesize/(2*3600) + 0.5*distance_between_source/3600.0]

        sample = []
        for i in range(grid):
            for j in range(grid):
                ra = start[0] + distance_between_source/3600. * i # degree
                dec = start[1] + distance_between_source/3600. * j

                # for easy transformation
                #c = coordinates.SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='fk5', equinox='J2000') 
                #pos = 'J2000 ' + str(c.to_string('hmsdms')) # unicode to str
                # Error: ancient astropy in ALMA Cluster!

                # flux scalling from PB
                theta_from_the_center = np.sqrt((ra - center[0])**2 + (dec - center[1])**2) # in degree
                scaled_flux = flux * np.exp(-np.log(2) * (2*theta_from_the_center/PB_in_deg)**2)

                pos = 'J2000 ' + str(ra) + 'deg ' + str(dec) + 'deg'

                sample.append([pos, scaled_flux, 'Jy', freq, 'point']) 
        
        # write to file
        with open(outputfile, 'w') as ofile:
            for item in sample:
                for i in item:
                    ofile.write("%s " % i)
                ofile.write('\n')

        return sample # ra, dec



    def read_sample(self, inputfile='sample.dat'):
        """
        Function to read the sample from ascii file
        if needed
        """

        sample = []
        with open(inputfile, 'r') as ifile:
            for line in ifile:
                item = line.strip().split()
                pos = item[0] + " " + item[1] + " " + item[2]
                flux = float(item[3])
                sample.append([pos, flux, item[4], item[5], item[6]])
                
        return sample



    def add_comp(self, ms_name, injected_source_list=[], 
        complist_name='injected_source.cl', resultms='injected.ms'):
        """ 
        split the original MS -- Data column
        make a component list
        add the fake source to the new ms
        split the MS again -- Corrected column
        """

        # remove complist if exist
        if os.path.exists(complist_name):
            os.system('rm -rf '+complist_name)

        copyms='injection.ms'
        if os.path.exists(copyms): # safety reason
            os.system('rm -rf '+copyms)
        # split the MS -- DATA COLUMN
        split(vis=ms_name,outputvis=copyms,keepmms=True,field="0",spw="",scan="",antenna="",correlation="",timerange="",
            intent="",array="",uvrange="",observation="",feed="",datacolumn="data",keepflags=True,width=1,timebin="0s",combine="")

        # add component
        cl.done()
        for i, pars in enumerate(injected_source_list):
            cl.addcomponent(dir=pars[0], flux=pars[1], fluxunit=pars[2], freq=pars[3], shape=pars[4])
        
        cl.rename(complist_name)
        cl.done()

        # fourier tranform the component
        ft(vis=copyms, complist=complist_name)

        # add the component
        uvsub(vis=copyms, reverse=True)

        # if resultms already exist
        if os.path.exists(resultms):
            os.system('rm -rf '+resultms)

        # split the MS -- CORRECTED COLUMN
        split(vis=copyms,outputvis=resultms,keepmms=True,field="0",spw="",scan="",antenna="",correlation="",timerange="",
            intent="",array="",uvrange="",observation="",feed="",datacolumn="corrected",keepflags=True,width=1,timebin="0s",combine="")

        # remove temporary MS
        os.system("rm -rf "+copyms)

        print "New MS: ", resultms

        return resultms
    

    def cleaning(self, ms_name='injected.ms', imgname='injected.ms.cont', niter=1000, threshold='', psfmode='clark', 
        interactive=False, mask='', imsize='', cell='', phasecenter='', weighting='briggs', robust=0.5, pbcor=False, regionrms=''):
        """ clean, calculate rms, export fits """

        # if image already exist
        for ifile in glob.glob(imgname+".*"):
            print "Removing previous image file:", ifile
            os.system("rm -rf "+ifile)

        print "CLEANing..."
        clean(vis=ms_name,imagename=imgname,outlierfile="",field="0",spw="",selectdata=True,
            timerange="",uvrange="",antenna="",scan="",observation="",intent="",mode="mfs",resmooth=False,gridmode="",
            wprojplanes=-1,facets=1,cfcache="cfcache.dir",rotpainc=5.0,painc=360.0,aterm=True,psterm=False,mterm=True,
            wbawp=False,conjbeams=True,epjtable="",interpolation="linear",niter=niter,gain=0.1,threshold=threshold,
            psfmode=psfmode,imagermode="csclean",ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=[],
            negcomponent=-1,smallscalebias=0.6,interactive=interactive,mask=mask,nchan=-1,start=0,width=1,outframe="",
            veltype="radio",imsize=imsize,cell=cell,phasecenter=phasecenter,restfreq="",stokes="I",weighting=weighting,
            robust=robust,uvtaper=False,outertaper=[''],innertaper=['1.0'],modelimage="",restoringbeam=[''],pbcor=pbcor,
            minpb=0.2,usescratch=False,noise="1.0Jy",npixels=0,npercycle=100,cyclefactor=1.5,cyclespeedup=-1,nterms=1,
            reffreq="",chaniter=False,flatnoise=True,allowchunk=False)

        # if fits file exists
        if os.path.exists(imgname+".fits"):
            os.system("rm rf "+imgname+".fits")

        exportfits(imagename=imgname+".image",fitsimage=imgname+".fits",velocity=False,
            optical=False,bitpix=-32,minpix=0,maxpix=-1,overwrite=False,dropstokes=False,stokeslast=True,history=True,
            dropdeg=False)


        print "Calculate RMS..."
        imstat(imagename=imgname+".image",axes=-1,region=regionrms,box="",chans="",stokes="I",listit=True,verbose=True,
            mask="",stretch=False,logfile=imgname+'.stat',append=True,algorithm="classic",fence=-1,center="mean",lside=True,zscore=-1,
            maxiter=-1,clmethod="auto")

        print "Result: ", imgname+".fits"

        return imgname



if __name__ == '__main__':
    ##
    # Run: 
    # casa -c completeness.py --nologger --log2term
    ##

    cT = completenessTest()
    
    #### a single value of flux ###
    # sample = cT.generate_sample(flux=1.02e-04, freq='241.23GHz', center=[65.8158364, -1.3425182], PB=24.0, beam=0.2, 
    #     imagesize=48.0, spacing=10, outputfile="sample.dat")
    
    # cT.add_comp("all_withweight_J0423-0120.ms", sample)
    
    # cT.cleaning(ms_name='injected.ms', imgname='injected_flux_3sigma.cont', niter=1000, threshold="0.035mJy", mask="circle[[65.8158364deg, -1.3425182deg], 5arcsec]", 
    #     imsize=1250, cell="0.04arcsec", phasecenter="J2000 04:23:15.800730 -01.20.33.065501", 
    #     regionrms="annulus[[65.8158364deg, -1.3425182deg], [5arcsec, 15arcsec]]")


    #### grid in flux ###
    fluxmin = 6.8e-05
    fluxmax = 3.4e-04
    nflux   = 16
    for flux in np.linspace(fluxmin, fluxmax, nflux):
        sample = cT.generate_sample(flux=flux, freq='241.23GHz', center=[65.8158364, -1.3425182], PB=24.0, beam=0.2, 
            imagesize=48.0, spacing=10, outputfile="sample.dat")
        
        cT.add_comp("all_withweight_J0423-0120.ms", sample)
        
        cT.cleaning(ms_name='injected.ms', imgname='injected_flux_'+ str(flux) +'.cont', niter=1000, threshold="0.035mJy", mask="circle[[65.8158364deg, -1.3425182deg], 5arcsec]", 
            imsize=1250, cell="0.04arcsec", phasecenter="J2000 04:23:15.800730 -01.20.33.065501", 
            regionrms="annulus[[65.8158364deg, -1.3425182deg], [5arcsec, 15arcsec]]")
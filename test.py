###########################################################
# directories

Dir = '/vol/alcina/data1/rfink/ngc3256_YMCs/'
scriptDir = Dir + 'flux/'
regionDir = Dir + 'flux/regions/'
imageDir_band3 = Dir + 'flux/'

sourcefile_band3 = regionDir + 'region_init'
source_regions = '/vol/alcina/data1/rfink/ngc3256_YMCs/flux/region_init_1'
imagefile_band3 = imageDir_band3 + 'member.uid___A001_X2d20_X3a88.NGC_3256_sci.spw5_17_19_23.cont.I.tt0.pbcor.fits'

outfile = scriptDir + 'region_imfit.crtf'
outfile_fitted = scriptDir + 'fittedregion_imfit.crtf'
###########################################################
# basic settings
import os, sys
ratio = 2.0


###########################################################
# function



def beam_get(imagename, regions, ratio=2.0, **kwarg):
    '''
    parameters
    imagename: The path to the CASA image file
    regions: The path to the CASA region file or the CASA
    region string.
    ratio: The ratio to be multiplied by the fitted ellipse to
    get the final elliptical aperture shape.
    **kwarg: other parameters that goes into the imfit
    '''
    beam=imfit(imagename=imagename,region=regions,model = 'gausfit_region.image',logfile='log_region_fit.txt',append = True,  **kwarg)
    regions = []
    for i in range(beam['results']['nelements']):
        component = 'component' + str(i)
        x_value=beam['results'][component]['shape']\
                 ['direction']['m0']['value']
        y_value=beam['results'][component]['shape']\
                 ['direction']['m1']['value']
        bmaj_value=beam['results'][component]\
                    ['shape']['majoraxis']['value']
        bmin_value=beam['results'][component]['shape']\
                    ['minoraxis']['value']
        pa_value=beam['results'][component]['shape']\
                  ['positionangle']['value']
        x=str(x_value)+'rad'
        y=str(y_value)+'rad'
        bmaj=str(bmaj_value/2.0*ratio)+'arcsec'
        bmin=str(bmin_value/2.0*ratio)+'arcsec'
        pa=str(pa_value)+'deg'
        region='ellipse[['+x+','+y+'],['+bmaj+','+bmin+'],'+pa+']'
        regions.append(region)

    return regions

def beam_to_regions(beam, ratio=2.0):
    '''
    Write the CASA beam object into CASA
    ------
    Parameters:
    beam: CASA beam object
    The beam acquired from task such as imfit
    ratio: float
    The ratio value to be multiplied to the beam object 
    '''
    regions = []
    for i in range(beam['results']['nelements']):
        component = 'component' +int(str(i))
        x_value=beam['results'][component]['shape']\
                 ['direction']['m0']['value']
        y_value=beam['results'][component]['shape']\
                 ['direction']['m1']['value']
        bmaj_value=beam['results'][component]\
                    ['shape']['majoraxis']['value']
        bmin_value=beam['results'][component]['shape']\
                    ['minoraxis']['value']
        pa_value=beam['results'][component]['shape']\
                  ['positionangle']['value']
        x=str(x_value)+'rad'
        y=str(y_value)+'rad'
        bmaj=str(bmaj_value/2.0*ratio)+'arcsec'
        bmin=str(bmin_value/2.0*ratio)+'arcsec'
        pa=str(pa_value)+'deg'
        region='ellipse[['+x+','+y+'],['+bmaj+','+bmin+'],'+pa+']'
        regions.append(region)

    return regions


###################################
#main  program

#beam_get(imagename='image_casa100ghz.image', regions='region_single', ratio=2.0)

lats_err = []; lons_err = []

regions = []
fittedRegions = []
fittedBeams = []

with open(source_regions, 'r') as infile:
    line = infile.readline()
    while line !='':
        line = infile.readline()
        regions.append(line)
regions.remove('')

with open (outfile, 'w') as out:
        out.write('#CRTFv0 CASA Region Text Format version 0\n')
        for line in regions:
            out.write(line)

for i, region in enumerate(regions):
    fittedBeam = beam_get(imagename='image_casa100ghz.image', regions=region, ratio=2.0)
    print (fittedBeam)
    fittedRegion = beam_to_regions(beam=fittedBeam, ratio = 2.0)
    fittedRegions += fittedRegion

with open (outfile_fitted, 'w') as outfile:
        outfile.write('#CRTFv0 CASA Region Text Format version 0\n')
        for line in fittedRegions:
            outfile.write(line)

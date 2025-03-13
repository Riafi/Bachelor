
###########################################################
# directories

Dir = '/vol/alcina/data1/rfink/ngc3256_YMCs/'
scriptDir = Dir + 'analyse/'
regionDir = Dir + '/analyse/regions/'
imageDir_band3 = Dir + '/analyse/cont_100GHz/'
imageDir_band6 = Dir + '/analyse/cont_345GHz/'
#imageDir_band3_2016 = Dir + '2016/Band3/'

sourcefile_band3 = regionDir + 'source_bandboth_3rms_noimfit.crtf'
imagefile_band3 = imageDir_band3 + 'ngc3256_band3_cont_12mext.pbcor'

sourcefile_band6 = regionDir + 'source_bandboth_3rms_noimfit.crtf'
imagefile_band6 = imageDir_band6 + 'ngc3256_band3_cont_345GHz_12mext.pbcor'

#sourcefile_band3_2016 = regionDir + 'starcluster_band3_2016.crtf'
#imagefile_band3_2016 = imageDir_band3_2016 + \
#                       'ngc40389overlap_band3_uvrange_robust_2_smooth.pbcor'

###########################################################
# basic settings

ratio = 2.0

bands = ['band 3', 'band 6']
continuum = dict.fromkeys(bands)

# band 3 image parameters for 2018 data
continuum['band 3'] = {}
continuum['band 3']['imagename'] = imagefile_band3
continuum['band 3']['region'] = sourcefile_band3
continuum['band 3']['outfile'] = regionDir + 'source_bandboth3_imfit_3rms.crtf'
continuum['band 3']['pbimage'] = imageDir_band3 + 'ngc3256_band3_cont_12mext.pb'
continuum['band 3']['rms'] = None

# band 6 image parameters for 2018 data.
continuum['band 6'] = {}
continuum['band 6']['imagename'] = imagefile_band6
continuum['band 6']['region'] = sourcefile_band6
continuum['band 6']['outfile'] = regionDir + 'source_bandboth6_imfit_3rms.crtf'
continuum['band 6']['pbimage'] = imageDir_band6 + 'ngc3256_band3_cont_345GHz_12mext.pb'
continuum['band 6']['rms'] = None

###########################################################
# function

def beam_get(imagename, region_init, ratio=2.0, **kwarg):
    '''
    parameters
    imagename: The path to the CASA image file
    regions: The path to the CASA region file or the CASA
    region string.
    ratio: The ratio to be multiplied by the fitted ellipse to
    get the final elliptical aperture shape.
    **kwarg: other parameters that goes into the imfit
    '''
    beam=imfit(imagename=imagename,region=region_init, **kwarg)
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

###########################################################
# main program 

# Set the default value of the ratio. 
for band in bands:
    continuum[band]['ratio'] = ratio
    if band in ['band3_p9_multiple', 'band3_p9_1b','HST_Pbeta','HST_Pbeta2']:
        continuum[band]['ratio'] = 2.0

# set the double gaussian fitting
for band in bands:
    continuum[band]['kwarg'] = {}
    if band == 'band 1':
        continuum[band]['kwarg']['estimates'] = scriptDir + 'estimate/band3_estimates.txt'


bands = ['band 3', 'band 6']
#bands = ['band7_2016_p5']
#bands = ['band3_2016_p5']
#bands = ['band3_2016_2_p201', 'band7_2016_2_p201']
#bands = ['HST_Pbeta']
#bands = ['band3_2016_p5', 'band7_2016_p5_band3Ap']
#bands = ['HST_Pbeta2']
#bands = ['HST_I']
#bands =['band7_2016_band3Ap']
#bands = ['HST_Pbeta_contsub']

lats_err = []; lons_err = []

for band in bands:
    regions = []
    with open(continuum[band]['region'], 'r') as infile:
        line = infile.readline()
        while line !='':
            line = infile.readline()
            regions.append(line)
    regions.remove('')

    # fit different regions with 2D Gaussian function and decide
    # the aperture size to measure the flux.
    fittedRegions = []
    fittedBeams = []
    for i, region in enumerate(regions):
#         rms = continuum[band]['rms'] / continuum[band]['pbcor'][i]
        fittedBeam = imfit(imagename = continuum[band]['imagename'],
                            region = region,
#                             rms = rms,
                            **continuum[band]['kwarg'])

        fittedRegion = beam_to_regions(fittedBeam, ratio=continuum[band]['ratio'])
        fittedRegions += fittedRegion
        if 'deconvolved' in fittedBeam and 'component0' in fittedBeam['deconvolved'] and 'shape' in fittedBeam['deconvolved']['component0'] and 'majoraxis' in fittedBeam['deconvolved']['component0']['shape']:
            print('deconvolved',',',fittedBeam['deconvolved']['component0']['shape']['majoraxis']['value'],',',fittedBeam['deconvolved']['component0']['shape']['majoraxiserror']['value'],',',fittedBeam['deconvolved']['component0']['shape']['minoraxis']['value'],',',fittedBeam['deconvolved']['component0']['shape']['minoraxiserror']['value'])
        else:

            if 'results' in fittedBeam and 'component0' in fittedBeam['results']:
                print('not deconvolved',',',fittedBeam['results']['component0']['shape']['majoraxis']['value'],',',fittedBeam['results']['component0']['shape']['majoraxiserror']['value'],',',fittedBeam['results']['component0']['shape']['minoraxis']['value'],',',fittedBeam['results']['component0']['shape']['minoraxiserror']['value'])
            else:
                print('Fit could not converge')
#         print(fittedBeam['deconvolved']['component0']['shape']['majoraxis'])
#         print(fittedBeam['deconvolved']['component0']['shape']['majoraxiserror'])
#         lats_err.append(fittedBeam['results']['component0']['shape']['direction']['error']['latitude']['value'])
#         lons_err.append(fittedBeam['results']['component0']['shape']['direction']['error']['longitude']['value'])

    # Export the regions into a file.
    with open (continuum[band]['outfile'], 'w') as outfile:
        outfile.write('#CRTFv0 CASA Region Text Format version 0\n')
        for line in fittedRegions:
            outfile.write(line+', color=green\n')


print(lats_err)
print(lons_err)

for band in bands:
    regions = []
    with open(continuum[band]['region'], 'r') as infile:
        line = infile.readline()
        while line !='':
            line = infile.readline()
            regions.append(line)
    regions.remove('')

    # fit different regions with 2D Gaussian function and decide
    # the aperture size to measure the flux.
    fittedRegions = []
    fittedBeams = []
    for i, region in enumerate(regions):
#         rms = continuum[band]['rms'] / continuum[band]['pbcor'][i]
        beam = imfit(imagename = continuum[band]['imagename'],
                            region = region,
#                             rms = rms,
                            **continuum[band]['kwarg'])

        #print the deconvolved parameter derived from imfit, if deconvolving fails, print the convolved result parameter derived from imfit
        if 'deconvolved' in beam and 'component0' in beam['deconvolved'] and 'peak' in beam['deconvolved']['component0']:
            print('deconvolved', ',',beam['deconvolved']['component0']['shape']['direction']['m0']['value'],',',beam['deconvolved']['component0']['shape']['direction']['m1']['value'],',',beam['deconvolved']['component0']['shape']['positionangle']['value'],',',beam['deconvolved']['component0']['shape']['direction']['error']['longitude']['value'],',',beam['deconvolved']['component0']['shape']['direction']['error']['latitude']['value'],',',beam['deconvolved']['component0']['flux']['value'][0],',',beam['deconvolved']['component0']['flux']['error'][0], ',',beam['deconvolved']['component0']['peak']['value'], ',' ,beam['deconvolved']['component0']['peak']['error'])
        else:

            if 'results' in beam and 'component0' in beam['results'] and 'peak' in beam['results']['component0']:
                print('not deconvolved',',', beam['results']['component0']['shape']['direction']['m0']['value'],',',beam['results']['component0']['shape']['direction']['m1']['value'],',',beam['results']['component0']['shape']['positionangle']['value'],',',beam['results']['component0']['shape']['direction']['error']['longitude']['value'],',',beam['results']['component0']['shape']['direction']['error']['latitude']['value'],',',beam['results']['component0']['flux']['value'][0],',',beam['results']['component0']['flux']['error'][0], ',',beam['results']['component0']['peak']['value'], ',' ,beam['results']['component0']['peak']['error'])
            else:
                print('Fit could not converge')

#print the deconvolved parameter derived from imfit, if deconvolving fails, print the convolved result parameter derived from imfit
print('for the estimate file')
for band in bands:
    regions = []
    with open(continuum[band]['region'], 'r') as infile:
        line = infile.readline()
        while line !='':
            line = infile.readline()
            regions.append(line)
    regions.remove('')

    # fit different regions with 2D Gaussian function and decide
    # the aperture size to measure the flux.
    fittedRegions = []
    fittedBeams = []
    for i, region in enumerate(regions):
#         rms = continuum[band]['rms'] / continuum[band]['pbcor'][i]
        beam = imfit(imagename = continuum[band]['imagename'],
                            region = region,
#                             rms = rms,
                            **continuum[band]['kwarg'])
        if 'deconvolved' in beam and 'component0' in beam['deconvolved'] and 'peak' in beam['deconvolved']['component0']:
            print(beam['deconvolved']['component0']['peak']['value'], ',' ,beam['deconvolved']['component0']['shape']['direction']['m0']['value'],',',beam['deconvolved']['component0']['shape']['direction']['m1']['value'],',',str(beam['deconvolved']['component0']['shape']['majoraxis']['value'])+'arcsec',',',str(beam['deconvolved']['component0']['shape']['minoraxis']['value'])+'arcsec',',',str(beam['deconvolved']['component0']['shape']['positionangle']['value'])+'deg')
        else:

            if 'results' in beam and 'component0' in beam['results'] and 'peak' in beam['results']['component0']:
                print(beam['results']['component0']['peak']['value'], ',' , beam['results']['component0']['shape']['direction']['m0']['value'],',',beam['results']['component0']['shape']['direction']['m1']['value'],',',str(beam['results']['component0']['shape']['majoraxis']['value'])+'arcsec',',',str(beam['results']['component0']['shape']['minoraxis']['value'])+'arcsec',',',str(beam['results']['component0']['shape']['positionangle']['value'])+'deg')
            else:
                print('Fit could not converge')

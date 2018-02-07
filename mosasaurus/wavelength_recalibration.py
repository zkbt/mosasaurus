# Zach, I am so sorry this code is so convoluted. See if you can make it work and if not let me know and I will try to help out. It may have become less compatible with your updates...
# To be run after reduce.py, but before show.py
# If you run mosasaurus in an interavtive session you can then copy and paste this whole code in there and then cross your fingers and hope
# If you want to read through in order of what gets executed, please scross down until you see "START HERE!"
# I tried to add comments to everything but I did this rather quickly so please let me know if something is not clear.

import sys
import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial import Legendre


def align_lines(stars, line_range, offset, line_pos, line_poses, plot=False):
    '''
    This function takes
        stars:    the extracted*.npy file
        line_range:    the pre-determined fixed alignment range for a given feature to use for the alignment
        offset:   how many angstroms to shift over when performing the cross-correlation of a given feature in a spectrum
        line_pos:    a list that this function will append values to; amount of shift in wavelength space for each star width
        line_poses:    a list of line_pos lists; amount of shifts in wavelength space for each star with all its widths
        plot:  do you want plots to pop up? (T/F)

    This code assumes that the stellar spectra may shift slowly throughout the night but do not suddenly jump around from one exposure to the next.

    '''

    # find the wavelengths between the wavelength range of the feature for the master star
    idxmaster = (starmaster['wavelength']>=line_range[0])*(starmaster['wavelength']<=line_range[1])
    # get the corresponding values from raw_counts; this will be the main thing to correlate against
    corrmaster = starmaster[width]['raw_counts'][idxmaster]
    # this is where the "true" value of the line is; i.e., the reference position on the master spectrum
    line = starmaster['wavelength'][idxmaster][0]
    # where that line falls in pixel space; uses the wavelength solution from mosasaurus
    linepx = np.interp(line, starmaster['wavelength'], starmaster['w'])

    # list of stars is the stars from "aperture" in mosasaurus, plus the master star appended at the end
    for s in range(len(stars)):

        if plot == True: print('star', s)

        # find the wavelengths between the wavelenth range of the feature for this particular star
        idxstar = (stars[s]['wavelength']>=line_range[0])*(stars[s]['wavelength']<=line_range[1])
        # this is the spectrum of the star that you want to correlate to the master
        corrstar = stars[s][width]['raw_counts']
        # need to know now many wavelengths are being covered
        arraylen = len(np.where(idxmaster)[0])
        # this is the pixel where we will start the shift
        initpx = np.where(idxstar)[0][0] - offset
        #if corrstar[initpx] == 0.0:
        #    initpx = np.where(corrstar != 0.)[0][0]
        #    corrmaster = wavemaster[initpx+offset:initpx+offset+arraylen]
        # empty lists for the correlation coeffients and the actual pixel shifts
        corrs, shifts = [], []

        # set the number of shifts you will try in the cross correlation; the finer the shift the more accurate everything will be but it will take longer
        for shift in np.linspace(0, arraylen+offset, (arraylen+offset)*10+1):
            # this is to make sure we're not losing counts at the ends of pixels when we're shifting by sub-pixels
            newbit = shift%1
            therest = 1 - newbit
            startpx = int(np.floor(shift))
            # create the array at the sub-pixel level that will be compared to the master star feature
            corrarray = corrstar[initpx+startpx : initpx+startpx+arraylen]*therest + corrstar[initpx+startpx+1 : initpx+startpx+arraylen+1]*newbit
            #if shift == 0:
            #    plt.plot(corrmaster)
            #    plt.plot(corrarray)
            #    plt.show()
            corrs.append(np.corrcoef(corrmaster, corrarray)[0,1])
            shifts.append(corrarray)

        if plot == True and (s != len(stars)-1):
            # can inspect where the code thinks the peak of the correlation is (careful! it doesn't always pick the right one!)
            print('1st corr')
            plt.plot(corrs)
            plt.axvline(np.where(corrs == np.max(corrs))[0])
            plt.title(str(line_range))
            plt.show()

        # try a shift based on the correlation coefficients
        # this first try may be wrong so we have to compare to past shifts
        firstshiftind = np.where(corrs == np.max(corrs))[0][0]
        firstrange = np.linspace(0, arraylen+offset, (arraylen+offset)*10+1)
        firstshift = firstrange[firstshiftind]
        # offset should be reasonably small since we assume the spectral features do not suddenly jump around from one exposure to the next
        if len(line_poses) < 13:
            while ((firstshift - offset) > 5.) or ((firstshift - offset) < -5.):
                print('extra corr', line_range, firstshift-offset)
                currentmaxind = np.where(corrs == np.max(corrs))
                # delete the correlation coefficient that is the maximum but is not providing a reasonably small shift
                corrs = np.delete(corrs, currentmaxind)
                # find a new shift value
                firstshiftind = np.where(corrs == np.max(corrs))[0][0]
                firstshift = firstrange[firstshiftind]
                if ((firstshift - offset) < 5.) and ((firstshift - offset) > -5.) and plot == True and s != len(stars)-1:
                    plt.plot(corrs)
                    plt.axvline(np.where(corrs == np.max(corrs))[0])
                    plt.title(str(line_range))
                    plt.show()

        # once we have made enough reasonably small shifts we can use past shifts to determine whether or not the next shift is a reasonable jump
        # there is probably a better way to do this... like going through the whole thing and then looking for outliers
        else:
            # look at the last 13 shifts in line position for this particular star
            med = np.median(np.array(line_poses)[:, s][-13:] - np.array(line_poses)[:, -1][-13:])
            std = np.std(np.array(line_poses)[:, s][-13:] - np.array(line_poses)[:, -1][-13:])
            # this is apparently a good number of sigma to go by
            timesstd = 22
            # eventually the value that gets save in line_poses is a wavelength; make the transformation into wavelength space so that we can compare the current shift to the previous ones
            if (pxtowavemaster(linepx + firstshift - offset)-line > (med + timesstd*std)) or (pxtowavemaster(linepx + firstshift - offset)-line < (med - timesstd*std)):
                print('range: ', med - timesstd*std, '-', med + timesstd*std, pxtowavemaster((linepx + firstshift - offset))-line)
                print('stddev corr', line_range)
                corrind = np.where((pxtowavemaster(linepx + firstrange - offset)-line >= (med - timesstd*std)) & (pxtowavemaster(linepx + firstrange - offset)-line <= (med + timesstd*std)))
                corrsclipped = np.array(corrs)[corrind]
                firstshiftind = np.where(corrsclipped == np.max(corrsclipped))[0][0] + corrind[0][0]
                firstshift = firstrange[firstshiftind]
                if plot == True and s != len(stars)-1:
                    plt.plot(corrs)
                    plt.axvline(np.where(corrsclipped == np.max(corrsclipped))[0][0] + corrind[0][0])
                    plt.title(str(line_range))
                    plt.show()

        # interpolate shift using a parabola to find the true peak cross correlation value
        parabinds = range(firstshiftind-2,firstshiftind+3)
        parabrange = np.linspace(parabinds[0], parabinds[-1], 100*(len(parabinds)-1))
        parabfit = np.polyfit(parabinds, np.array(corrs)[parabinds], 2)
        parabfit1d = np.poly1d(parabfit)
        parabcorrs = parabfit1d(parabrange)
        parabshiftind = np.where(parabcorrs == np.max(parabcorrs))[0]
        parabshift = parabrange[parabshiftind]
        #if plot == True:
        #    plt.plot(parabinds, np.array(corrs)[parabinds])
        #    plt.plot(parabrange, parabcorrs)
        #    plt.axvline(parabshift)
        #   plt.show()
        frac = parabshift - firstshiftind
        finalshift = firstshift + frac*np.median(np.diff(firstrange))

        # transform all these shifts in pixel space into wavelength space
        newlinepx = (linepx + finalshift - offset)[0]
        newline = pxtowavemaster(newlinepx)
        #print newline
        # remember that the last "star" in the list is actually just the master spectrum; again, not sure why I did it this way
        if s == len(stars)-1: line_pos.append(line)
        else: line_pos.append(newline)
        if plot == True: print(finalshift - offset)
        #print newline

        # this was left-over from some test I did to see how much the GJ1132 spectrum was stretching/shifting throughout the night
        #shift1132.append((starline+finalshift-offset)[0])

        #newbit = finalshift%1
        #therest = 1 - newbit
        #startpx = int(np.floor(finalshift))
        #newarray = corrstar[initpx+startpx:initpx+startpx+arraylen]*therest + corrstar[initpx+startpx+1:initpx+startpx+arraylen+1]*newbit

        #if plot == True:
        #    plt.plot(newarray, label=str(s))

    #if plot == True:
    #    plt.plot(corrmaster, label='corrmaster')
    #    plt.title(str(line_range))
    #    plt.legend()
    #    plt.show()

################ START HERE!! ######################

# 2016_02_27-28
mastern = '1825'
starmasterstr = '/h/mulan0/data/working/GJ1132b_ut160227_28/multipleapertures/aperture_544_1024/extracted'+mastern+'.npy'
# get the values to recreate the wavelength solution from mosasaurus
coef, domain =  np.load('/h/mulan0/data/working/GJ1132b_ut160227_28/multipleapertures/aperture_544_1024/aperture_544_1024_wavelengthcalibration.npy')[()]
# 2016_03_03-04
#mastern = '0945'
#starmaster = np.load('/h/mulan0/data/working/GJ1132b_ut160303_04/multipleapertures/aperture_588_1045/extracted'+mastern+'.npy')[()]
#coef, domain =  np.load('/h/mulan0/data/working/GJ1132b_ut160303_04/multipleapertures/aperture_588_1045/aperture_588_1045_wavelengthcalibration.npy')[()]
# 2016_03_08-09
#mastern = '0348'
#starmasterstr = '/h/mulan0/data/working/GJ1132b_ut160308_09/multipleapertures/aperture_587_1042/extracted'+mastern+'.npy'
#coef, domain =  np.load('/h/mulan0/data/working/GJ1132b_ut160308_09/multipleapertures/aperture_587_1042/aperture_587_1042_wavelengthcalibration.npy')[()]
# 2016_04_16-17
#mastern = '1252'
#starmasterstr = '/h/mulan0/data/working/GJ1132b_ut160416_17/multipleapertures/aperture_584_1052/extracted'+mastern+'.npy'
#coef, domain =  np.load('/h/mulan0/data/working/GJ1132b_ut160416_17/multipleapertures/aperture_584_1052/aperture_584_1052_wavelengthcalibration.npy')[()]
# 2016_04_21-22
#mastern = '1386'
#starmaster = np.load('/h/mulan0/data/working/GJ1132b_ut160421_22/multipleapertures/aperture_587_1049/extracted'+mastern+'.npy')[()]
#oef, domain =  np.load('/h/mulan0/data/working/GJ1132b_ut160421_22/multipleapertures/aperture_587_1049/aperture_587_1049_wavelengthcalibration.npy')[()]

# this is a check to make sure you have commented in the correct star (this is obviously a hack that could be done better...)
response = raw_input('Is this the correct starmaster for this dataset? [y/n] \n     ' + starmasterstr + '\n')
if response == 'n': sys.exit('Please edit wavelength_recalibration.py to point to the correct starmaster')
elif response == 'y': starmaster = np.load(starmasterstr)[()]

# reacreate the wavelength solution (way to go from pixel space to wavelength space) from mosasaurus
pxtowavemaster = Legendre(coef, domain)
# apertures form a given night (Need to have run mosasaurus in ipython or similar and then copy and paste this monster code in. This is not ideal.)
apertures = r.mask.apertures
makeplot = True
UV_poses, O2_poses, Ca1_poses, Ca2_poses, Ca3_poses, H2O_poses = [], [], [], [], [], []
x_poses = []
#shifts_1132 = []
# Fixed alignment ranges for each prominent freature
align_UV = (6870, 6900)
align_O2 = (7580, 7650)
#align_H2Osmall = (8220, 8260)
align_Ca1 = (8490, 8525)
align_Ca2 = (8535, 8580)
align_Ca3 = (8650, 8700)
align_H2O = (9300, 9700)
for n in r.obs.nScience:
    extracted = '{0:04}.npy'.format(n)
    print('working on frame {0:04}'.format(n))
    stars = []
    for a in apertures:
        extractedpathname = a.directory + 'extracted' + extracted
        stars.append(np.load(extractedpathname)[()])
    # append to the list of stars the master spectrum; not sure why I did it this way but we do need to know the "master" line position in wavelength space so we can later calculate how many wavelengths to shift over
    stars.append(starmaster)

    # need to do each width separately
    for width in a.widths:

        # the change in position of each of these features
        UV_pos, O2_pos, Ca1_pos, Ca2_pos, Ca3_pos, H2O_pos = [], [], [], [], [], []
        # NOW GO UP TO THE aling_star FUNCTION
        #align_lines(stars, align_UV, 10, UV_pos, UV_poses, makeplot)
        align_lines(stars, align_O2, 20, O2_pos, O2_poses, makeplot)
        #H2Osmall_pos = []
        #align_lines(stars, align_H2Osmall, 20, H2Osmall_pos, makeplot)
        align_lines(stars, align_Ca1, 10, Ca1_pos, Ca1_poses, makeplot)
        align_lines(stars, align_Ca2, 10, Ca2_pos, Ca2_poses, makeplot)
        align_lines(stars, align_Ca3, 10, Ca3_pos, Ca3_poses, makeplot)
        align_lines(stars, align_H2O, 10, H2O_pos, H2O_poses, makeplot)

        #UV_poses.append(UV_pos)
        O2_poses.append(O2_pos)
        Ca1_poses.append(Ca1_pos)
        Ca2_poses.append(Ca2_pos)
        Ca3_poses.append(Ca3_pos)
        H2O_poses.append(H2O_pos)

        # do some list re-arranging; we need to know what the shift and stretch in wavelength are for each star
        # x_pos just means any of the feature positions we've calculated
        x_pos = np.array([O2_pos, Ca1_pos, Ca2_pos, Ca3_pos, H2O_pos])
        # here is where we remembered the master line position for each feature so we can subtract it off and know how much we moved
        x_poses.append((x_pos.transpose() - x_pos[:, -1]).transpose())
        known_wave = x_pos[:,-1]

        for i in range(len(stars)-1):
            # create a new list wavelength solution that goes from the wavelength mosasurus gives to the new wavelength we have moved to
            fit = np.polyfit(x_pos[:,i], known_wave, 1)
            fit1d = np.poly1d(fit)
            # save the new wavelength array, as well as the coefficients so that we can re-create this fine-tuned wavelength solution later on in detrendersaurus
            newwavelength = stars[i]['wavelength']
            stars[i][width]['wavelength_adjusted'] = fit1d(newwavelength)
            stars[i][width]['stretch'] = fit[0]
            stars[i][width]['shift'] = fit[1]
            a = apertures[i]
            extractedpathname = a.directory + 'extracted' + extracted
            np.save(extractedpathname, stars[i])

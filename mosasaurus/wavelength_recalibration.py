import sys
import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial import Legendre
def align_lines(stars, line_range, offset, line_pos, line_poses, plot=False):
    idxmaster = (starmaster['wavelength']>=line_range[0])*(starmaster['wavelength']<=line_range[1])
    corrmaster = starmaster[width]['raw_counts'][idxmaster]
    line = starmaster['wavelength'][idxmaster][0]
    #print line
    linepx = np.interp(line, starmaster['wavelength'], starmaster['w'])

    for s in range(len(stars)):

        if plot == True: print 'star', s

        idxstar = (stars[s]['wavelength']>=line_range[0])*(stars[s]['wavelength']<=line_range[1])
        corrstar = stars[s][width]['raw_counts']
        arraylen = len(np.where(idxmaster)[0])
        initpx = np.where(idxstar)[0][0] - offset
        #if corrstar[initpx] == 0.0: 
        #    initpx = np.where(corrstar != 0.)[0][0]
        #    corrmaster = wavemaster[initpx+offset:initpx+offset+arraylen]
        corrs = []
        shifts = []

        for shift in np.linspace(0, arraylen+offset, (arraylen+offset)*10+1):
            newbit = shift%1
            therest = 1 - newbit
            startpx = int(np.floor(shift))
            corrarray = corrstar[initpx+startpx:initpx+startpx+arraylen]*therest + corrstar[initpx+startpx+1:initpx+startpx+arraylen+1]*newbit
            #if shift == 0:
            #    plt.plot(corrmaster)
            #    plt.plot(corrarray)
            #    plt.show()
            corrs.append(np.corrcoef(corrmaster, corrarray)[0,1])
            shifts.append(corrarray)
        
        if plot == True and (s != len(stars)-1):
            print '1st corr'
            plt.plot(corrs)
            plt.axvline(np.where(corrs == np.max(corrs))[0])
            plt.title(str(line_range))
            plt.show()
        
        firstshiftind = np.where(corrs == np.max(corrs))[0][0]
        firstrange = np.linspace(0, arraylen+offset, (arraylen+offset)*10+1)
        firstshift = firstrange[firstshiftind]
        # either use need to make offsets small (if under 10 points) or use med and std of past shifts to constrain new shift
        if len(line_poses) < 13:
            while ((firstshift - offset) > 5.) or ((firstshift - offset) < -5.):
                print 'extra corr', line_range, firstshift-offset
                currentmaxind = np.where(corrs == np.max(corrs))
                corrs = np.delete(corrs, currentmaxind)
                firstshiftind = np.where(corrs == np.max(corrs))[0][0]
                firstshift = firstrange[firstshiftind]
                if ((firstshift - offset) < 5.) and ((firstshift - offset) > -5.) and plot == True and s != len(stars)-1:
                    plt.plot(corrs)
                    plt.axvline(np.where(corrs == np.max(corrs))[0])
                    plt.title(str(line_range))
                    plt.show()
    
        else: 
            #if len(line_poses) > 10:
            med = np.median(np.array(line_poses)[:, s][-13:] - np.array(line_poses)[:, -1][-13:])
            std = np.std(np.array(line_poses)[:, s][-13:] - np.array(line_poses)[:, -1][-13:])
            timesstd = 22
            if (pxtowavemaster(linepx + firstshift - offset)-line > (med + timesstd*std)) or (pxtowavemaster(linepx + firstshift - offset)-line < (med - timesstd*std)):
                print 'range: ', med - timesstd*std, '-', med + timesstd*std, pxtowavemaster((linepx + firstshift - offset))-line
                print 'stddev corr', line_range
                corrind = np.where((pxtowavemaster(linepx + firstrange - offset)-line >= (med - timesstd*std)) & (pxtowavemaster(linepx + firstrange - offset)-line <= (med + timesstd*std)))
                corrsclipped = np.array(corrs)[corrind]
                firstshiftind = np.where(corrsclipped == np.max(corrsclipped))[0][0] + corrind[0][0]
                firstshift = firstrange[firstshiftind]
                if plot == True and s != len(stars)-1:
                    plt.plot(corrs)
                    plt.axvline(np.where(corrsclipped == np.max(corrsclipped))[0][0] + corrind[0][0])
                    plt.title(str(line_range))
                    plt.show()
           
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
        

        newlinepx = (linepx + finalshift - offset)[0]
        newline = pxtowavemaster(newlinepx)
        #print newline
        if s == len(stars)-1: line_pos.append(line)
        else: line_pos.append(newline)
        if plot == True: print finalshift - offset
        #print newline

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


# to be run after reduce.py, but before show.py
# 2016_02_27-28
mastern = '1825'
starmasterstr = '/h/mulan0/data/working/GJ1132b_ut160227_28/multipleapertures/aperture_544_1024/extracted'+mastern+'.npy'
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

response = raw_input('Is this the correct starmaster for this dataset? [y/n] \n     ' + starmasterstr + '\n')
if response == 'n': sys.exit('Please edit wavelength_recalibration.py to point to the correct starmaster')
elif response == 'y': starmaster = np.load(starmasterstr)[()]

pxtowavemaster = Legendre(coef, domain)
apertures = r.mask.apertures
makeplot = False
UV_poses, O2_poses, Ca1_poses, Ca2_poses, Ca3_poses, H2O_poses = [], [], [], [], [], []
x_poses = []
#shifts_1132 = []
for n in r.obs.nScience:
    extracted = '{0:04}.npy'.format(n)
    print 'working on frame {0:04}'.format(n)
    stars = []
    for a in apertures:
        extractedpathname = a.directory + 'extracted' + extracted
        stars.append(np.load(extractedpathname)[()])
    stars.append(starmaster)
    for width in a.widths:
        #width = 6.0
        align_UV = (6870, 6900)
        align_O2 = (7580, 7650)
        #align_H2Osmall = (8220, 8260)
        align_Ca1 = (8490, 8525)
        align_Ca2 = (8535, 8580)
        align_Ca3 = (8650, 8700)
        align_H2O = (9300, 9700)

        UV_pos, O2_pos, Ca1_pos, Ca2_pos, Ca3_pos, H2O_pos = [], [], [], [], [], []
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

        x_pos = np.array([O2_pos, Ca1_pos, Ca2_pos, Ca3_pos, H2O_pos])
        x_poses.append((x_pos.transpose() - x_pos[:, -1]).transpose())
        known_wave = x_pos[:,-1]

        for i in range(len(stars)-1):
            fit = np.polyfit(x_pos[:,i], known_wave, 1)
            fit1d = np.poly1d(fit)
            newwavelength = stars[i]['wavelength']
            stars[i][width]['wavelength_adjusted'] = fit1d(newwavelength)
            stars[i][width]['stretch'] = fit[0]
            stars[i][width]['shift'] = fit[1]
            a = apertures[i]
            extractedpathname = a.directory + 'extracted' + extracted
            np.save(extractedpathname, stars[i])

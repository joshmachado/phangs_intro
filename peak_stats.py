from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import pandas as pd
import os.path
import scipy.stats


def retrieve(noise, source=[], res=[]):

        loc = '/Users/josh/projects/intro'
        if noise=='homogenized':
            for z in range(len(source)):
                mean_dist = np.zeros([4,3])
                mean_beam_sep = np.zeros([4,3])
                percentiles = np.zeros([4,15]) #5th, 16th, 50th, 84th, 95th
                perc = [5,16,50,84,95]
                for i in range(len(res)):
                    fp = '/Users/josh/projects/intro/'+str(source[z])+'/'+str(source[z])+'_12m+7m+tp_co21_'+str(res[i])+'pc_props.fits.bz2'
                    if os.path.isfile(fp):
                        tab = Table.read(fp)
                        cloudnum = np.array(tab['CLOUDNUM'])
                        distance = np.array(tab['DISTANCE_PC'])
                        x = np.array(tab['XCTR_DEG'])
                        y = np.array(tab['YCTR_DEG'])
                        xs = (x - np.mean(x))*np.cos(np.deg2rad(y)) #correction RA = RAcos(Dec)
                        dist = np.zeros((len(x),len(x)))
                        nn = np.zeros(len(x)) #stores CloudNum of nearest neighbor
                        mindist = np.zeros(len(x)) #stores distance to nearest neighbor
                        nn2 = np.zeros(len(x)) #stores CloudNum of 2nd nearest neighbor
                        mindist2 = np.zeros(len(x)) #stores distance to 2nd nearest neighbor
                        nn3 = np.zeros(len(x)) #stores CloudNum of 3rd nearest neighbor
                        mindist3 = np.zeros(len(x)) #stores distance to 3rd nearest neighbor
                        k = 0
                        j = 0
                        while k < len(x):
                            for j in range(len(x)):
                                dist[k,j] = np.sqrt(np.square(xs[k] - xs[j]) + np.square(y[k] - y[j]))
                            k+=1
                        dist[dist == 0] = np.nan
                        for k in range(len(x)):
                            ind = np.where(dist[k] == np.nanmin(dist[k])) #index of nearest neighbor
                            ind2 = np.where(dist[k] == np.partition(dist[k], 1)[1]) #index of 2nd nearest neighbor
                            ind3 = np.where(dist[k] == np.partition(dist[k], 2)[2]) #index of 3rd nearest neighbor
                            mindist[k] = np.deg2rad(dist[k,int(ind[0])])*distance[0]
                            nn[k] = int(ind[0])+1
                            mindist2[k] = np.deg2rad(dist[k,int(ind2[0])])*distance[0]
                            nn2[k] = int(ind2[0])+1
                            mindist3[k] = np.deg2rad(dist[k,int(ind3[0])])*distance[0]
                            nn3[k] = int(ind3[0])+1


                        ### LETS GATHER SOME STATS ###
                        mean_dist[i] = [np.mean(mindist), np.mean(mindist2), np.mean(mindist3)]
                        mean_cloud_sep = [mindist, mindist2, mindist3]/np.nanmean(tab['RAD_PC'])
                        beam_sep = [mindist, mindist2, mindist3]/tab['BEAMFWHM_PC'][0]
                        percentiles[i] = [np.percentile(mindist, 5), np.percentile(mindist2, 5), np.percentile(mindist3, 5),
                        np.percentile(mindist, 16), np.percentile(mindist2, 16), np.percentile(mindist3, 16),
                        np.percentile(mindist, 50), np.percentile(mindist2, 50), np.percentile(mindist3, 50),
                        np.percentile(mindist, 84), np.percentile(mindist2, 84), np.percentile(mindist3, 84),
                        np.percentile(mindist, 95), np.percentile(mindist2, 95), np.percentile(mindist3, 95)] ##Go back and make this a loop
                        mean_beam_sep[i] = [np.mean(mindist), np.mean(mindist2), np.mean(mindist3)]/tab['BEAMFWHM_PC'][0]


                        ### Compile info into dataframe ###
                        cat = pd.DataFrame({'cloudnum':cloudnum, 'x':x, 'corr_x':xs, 'y':y,
                        'nn_index':nn,  'nn2_index':nn2, 'nn3_index':nn3,
                        'min_dist':mindist, 'min_dist2nd':mindist2, 'min_dist3rd':mindist3,
                        'mean_cloud_sep_nn':mean_cloud_sep[0], 'mean_cloud_sep_nn2':mean_cloud_sep[1], 'mean_cloud_sep_nn3':mean_cloud_sep[2],
                        'beam_sep_nn':beam_sep[0], 'beam_sep_nn2':beam_sep[1], 'beam_sep_nn3':beam_sep[2]})
                        cat.to_csv(loc+'/'+source[z]+'/'+source[z]+'_'+str(res[i])+'pc_cloud_stats.csv', index=False)
                        print(str(res[i])+'pc done')
                    print('Gathering stats...')

                    peak_stats = pd.DataFrame({'res_pc':res,
                    'mean_dist_nn':mean_dist[:,0], 'mean_dist_nn2':mean_dist[:,1], 'mean_dist_nn3':mean_dist[:,2],
                    '5th_nn':percentiles[:,0], '5th_nn2':percentiles[:,1], '5th_nn3':percentiles[:,2],
                    '16th_nn':percentiles[:,3], '16th_nn2':percentiles[:,4], '16th_nn3':percentiles[:,5],
                    '50th_nn':percentiles[:,6], '50th_nn2':percentiles[:,7], '50th_nn3':percentiles[:,8],
                    '84th_nn':percentiles[:,9], '84th_nn2':percentiles[:,10], '84th_nn3':percentiles[:,11],
                    '95th_nn':percentiles[:,12], '95th_nn2':percentiles[:,13], '95th_nn3':percentiles[:,14],
                    'mean_beam_sep':mean_beam_sep[:,0],'mean_beam_sep2':mean_beam_sep[:,1], 'mean_beam_sep3':mean_beam_sep[:,2]})

                    peak_stats.to_csv(loc+'/'+source[z]+'/'+source[z]+'_stats.csv', index=False)
        if noise=='matched':
            for z in range(len(source)):
                mean_dist = np.zeros([4,3])
                mean_beam_sep = np.zeros([4,3])
                percentiles = np.zeros([4,15]) #5th, 16th, 50th, 84th, 95th
                perc = [5,16,50,84,95]
                for i in range(len(res)):
                    fp = '/Users/josh/projects/intro/'+str(source[z])+'/matched/'+str(source[z])+'_12m+7m+tp_co21_'+str(res[i])+'pc_props.fits.bz2'
                    if os.path.isfile(fp):
                        tab = Table.read(fp)
                        cloudnum = np.array(tab['CLOUDNUM'])
                        distance = np.array(tab['DISTANCE_PC'])
                        x = np.array(tab['XCTR_DEG'])
                        y = np.array(tab['YCTR_DEG'])
                        xs = (x - np.mean(x))*np.cos(np.deg2rad(y)) #correction RA = RAcos(Dec)
                        dist = np.zeros((len(x),len(x)))
                        nn = np.zeros(len(x)) #stores CloudNum of nearest neighbor
                        mindist = np.zeros(len(x)) #stores distance to nearest neighbor
                        nn2 = np.zeros(len(x)) #stores CloudNum of 2nd nearest neighbor
                        mindist2 = np.zeros(len(x)) #stores distance to 2nd nearest neighbor
                        nn3 = np.zeros(len(x)) #stores CloudNum of 3rd nearest neighbor
                        mindist3 = np.zeros(len(x)) #stores distance to 3rd nearest neighbor
                        k = 0
                        j = 0
                        while k < len(x):
                            for j in range(len(x)):
                                dist[k,j] = np.sqrt(np.square(xs[k] - xs[j]) + np.square(y[k] - y[j]))
                            k+=1
                        dist[dist == 0] = np.nan
                        for k in range(len(x)):
                            ind = np.where(dist[k] == np.nanmin(dist[k])) #index of nearest neighbor
                            ind2 = np.where(dist[k] == np.partition(dist[k], 1)[1]) #index of 2nd nearest neighbor
                            ind3 = np.where(dist[k] == np.partition(dist[k], 2)[2]) #index of 3rd nearest neighbor
                            mindist[k] = np.deg2rad(dist[k,int(ind[0])])*distance[0]
                            nn[k] = int(ind[0])+1
                            mindist2[k] = np.deg2rad(dist[k,int(ind2[0])])*distance[0]
                            nn2[k] = int(ind2[0])+1
                            mindist3[k] = np.deg2rad(dist[k,int(ind3[0])])*distance[0]
                            nn3[k] = int(ind3[0])+1


                        ### LETS GATHER SOME STATS ###
                        mean_dist[i] = [np.mean(mindist), np.mean(mindist2), np.mean(mindist3)]
                        mean_cloud_sep = [mindist, mindist2, mindist3]/np.nanmean(tab['RAD_PC'])
                        beam_sep = [mindist, mindist2, mindist3]/tab['BEAMFWHM_PC'][0]
                        percentiles[i] = [np.percentile(mindist, 5), np.percentile(mindist2, 5), np.percentile(mindist3, 5),
                        np.percentile(mindist, 16), np.percentile(mindist2, 16), np.percentile(mindist3, 16),
                        np.percentile(mindist, 50), np.percentile(mindist2, 50), np.percentile(mindist3, 50),
                        np.percentile(mindist, 84), np.percentile(mindist2, 84), np.percentile(mindist3, 84),
                        np.percentile(mindist, 95), np.percentile(mindist2, 95), np.percentile(mindist3, 95)] ##Go back and make this a loop
                        mean_beam_sep[i] = [np.mean(mindist), np.mean(mindist2), np.mean(mindist3)]/tab['BEAMFWHM_PC'][0]


                        ### Compile info into dataframe ###
                        cat = pd.DataFrame({'cloudnum':cloudnum, 'x':x, 'corr_x':xs, 'y':y,
                        'nn_index':nn,  'nn2_index':nn2, 'nn3_index':nn3,
                        'min_dist':mindist, 'min_dist2nd':mindist2, 'min_dist3rd':mindist3,
                        'mean_cloud_sep_nn':mean_cloud_sep[0], 'mean_cloud_sep_nn2':mean_cloud_sep[1], 'mean_cloud_sep_nn3':mean_cloud_sep[2],
                        'beam_sep_nn':beam_sep[0], 'beam_sep_nn2':beam_sep[1], 'beam_sep_nn3':beam_sep[2]})
                        cat.to_csv(loc+'/'+source[z]+'/matched/'+source[z]+'_'+str(res[i])+'pc_cloud_stats.csv', index=False)
                        print(str(res[i])+'pc done')
                    print('Gathering stats...')

                    peak_stats = pd.DataFrame({'res_pc':res,
                    'mean_dist_nn':mean_dist[:,0], 'mean_dist_nn2':mean_dist[:,1], 'mean_dist_nn3':mean_dist[:,2],
                    '5th_nn':percentiles[:,0], '5th_nn2':percentiles[:,1], '5th_nn3':percentiles[:,2],
                    '16th_nn':percentiles[:,3], '16th_nn2':percentiles[:,4], '16th_nn3':percentiles[:,5],
                    '50th_nn':percentiles[:,6], '50th_nn2':percentiles[:,7], '50th_nn3':percentiles[:,8],
                    '84th_nn':percentiles[:,9], '84th_nn2':percentiles[:,10], '84th_nn3':percentiles[:,11],
                    '95th_nn':percentiles[:,12], '95th_nn2':percentiles[:,13], '95th_nn3':percentiles[:,14],
                    'mean_beam_sep':mean_beam_sep[:,0],'mean_beam_sep2':mean_beam_sep[:,1], 'mean_beam_sep3':mean_beam_sep[:,2]})

                    peak_stats.to_csv(loc+'/'+source[z]+'/matched/'+source[z]+'_stats.csv', index=False)



def pull_props(noise, source=[], res=[]):
    if noise=='homogenized':
        for j in range(len(res)):
            for i in range(len(source)):
                fp = '/Users/josh/projects/intro/'+str(source[i])+'/'+str(source[i])+'_12m+7m+tp_co21_'+str(res[j])+'pc_props.fits.bz2'
                stats = '/Users/josh/projects/intro/'+str(source[i])+'/'+str(source[i])+'_'+str(res[j])+'pc_cloud_stats.csv'
                if os.path.isfile(fp):
                    tab = Table.read(fp)
                    cat = pd.read_csv(stats)
                    df = pd.DataFrame(cat)
                    prop = ['MLUM_MSUN', 'SIGV_KMS', 'RAD_PC']
                    prop_nn = np.zeros(len(cat))
                    prop_nn2 = np.zeros(len(cat))
                    prop_nn3 = np.zeros(len(cat))
                    for z in range(len(prop)):
                        for k in range(len(cat)):
                            cloudnum = int(cat['cloudnum'][k])
                            nn = int(cat['nn_index'][cloudnum-1])
                            nn2 = int(cat['nn2_index'][cloudnum-1])
                            nn3 = int(cat['nn3_index'][cloudnum-1])
                            prop_nn[k] = tab[prop[z]][nn-1]
                            prop_nn2[k] = tab[prop[z]][nn2-1]
                            prop_nn3[k] = tab[prop[z]][nn3-1]
                        df[str(prop[z])+'_nn'] = prop_nn
                        df[str(prop[z])+'_nn2'] = prop_nn2
                        df[str(prop[z])+'_nn3'] = prop_nn3
                        df.to_csv('/Users/josh/projects/intro/'+str(source[i])+'/'+str(source[i])+'_'+str(res[j])+'pc_cloud_stats.csv', index=False)
    if noise=='matched':
        for j in range(len(res)):
            for i in range(len(source)):
                fp = '/Users/josh/projects/intro/'+str(source[i])+'/matched/'+str(source[i])+'_12m+7m+tp_co21_'+str(res[j])+'pc_props.fits.bz2'
                stats = '/Users/josh/projects/intro/'+str(source[i])+'/matched/'+str(source[i])+'_'+str(res[j])+'pc_cloud_stats.csv'
                if os.path.isfile(fp):
                    tab = Table.read(fp)
                    cat = pd.read_csv(stats)
                    df = pd.DataFrame(cat)
                    prop = ['MLUM_MSUN', 'SIGV_KMS', 'RAD_NODC']
                    prop_nn = np.zeros(len(cat))
                    prop_nn2 = np.zeros(len(cat))
                    prop_nn3 = np.zeros(len(cat))
                    for z in range(len(prop)):
                        for k in range(len(cat)):
                            cloudnum = int(cat['cloudnum'][k])
                            nn = int(cat['nn_index'][cloudnum-1])
                            nn2 = int(cat['nn2_index'][cloudnum-1])
                            nn3 = int(cat['nn3_index'][cloudnum-1])
                            prop_nn[k] = tab[prop[z]][nn-1]
                            prop_nn2[k] = tab[prop[z]][nn2-1]
                            prop_nn3[k] = tab[prop[z]][nn3-1]
                        df[str(prop[z])+'_nn'] = prop_nn
                        df[str(prop[z])+'_nn2'] = prop_nn2
                        df[str(prop[z])+'_nn3'] = prop_nn3
                        df.to_csv('/Users/josh/projects/intro/'+str(source[i])+'/matched/'+str(source[i])+'_'+str(res[j])+'pc_cloud_stats.csv', index=False)



def pull_pairs(source, res, noise, prop):
    loc = '/Users/josh/projects/intro/'
    global first_pairs, second_pairs, third_pairs, all_pairs
    props_dict = {
    "distance": ['min_dist', 'min_dist2nd', 'min_dist3rd'],
    "mean_cloud_sep": ['mean_cloud_sep_nn', 'mean_cloud_sep_nn2', 'mean_cloud_sep_nn3'],
    "beam_sep": ['beam_sep_nn', 'beam_sep_nn2', 'beam_sep_nn3'],
    "MLUM_MSUN": ['MLUM_MSUN_nn', 'MLUM_MSUN_nn2', 'MLUM_MSUN_nn3'],
    "SIGV_KMS": ['SIGV_KMS_nn', 'SIGV_KMS_nn2', 'SIGV_KMS_nn3'],
    "RAD_PC": ['RAD_PC_nn', 'RAD_PC_nn2', 'RAD_PC_nn3']
    }


    if noise=='homogenized':
        fp = loc+str(source)+'/'+str(source)+'_'+str(res)+'pc_cloud_stats.csv'
        if os.path.isfile(fp):
            cat = pd.read_csv(fp)
            tab = Table.read(loc+str(source)+'/'+str(source)+'_12m+7m+tp_co21_'+str(res)+'pc_props.fits.bz2')
            cloud = cat['cloudnum']
            nn_index, nn2_index, nn3_index = cat['nn_index'], cat['nn2_index'], cat['nn3_index']
            if prop == 'MLUM_MSUN' or prop=='SIGV_KMS' or prop=='RAD_PC':
                cloud_prop = tab[prop]
                nn_prop, nn2_prop, nn3_prop = cat[str(props_dict[prop][0])], cat[str(props_dict[prop][1])], cat[str(props_dict[prop][2])]
                first_pairs = pd.DataFrame({'cloud_index':cloud-1, 'nn_index':nn_index, 'cloud_prop':cloud_prop, 'nn_'+str(prop):nn_prop})
                second_pairs = pd.DataFrame({'cloud_index':cloud-1, 'nn2_index':nn2_index, 'cloud_prop':cloud_prop, 'nn2_'+str(prop):nn2_prop})
                third_pairs = pd.DataFrame({'cloud_index':cloud-1, 'nn3_index':nn3_index, 'cloud_prop':cloud_prop, 'nn3_'+str(prop):nn3_prop})
                all_pairs = pd.DataFrame({'cloud_index':cloud-1, 'nn_index':nn_index, 'nn2_index':nn2_index, 'nn3_index':nn3_index,
                'cloud_prop':cloud_prop, 'nn_'+str(prop):nn_prop, 'nn2_'+str(prop):nn2_prop, 'nn3_'+str(prop):nn3_prop})
                return first_pairs.dropna(), second_pairs.dropna(), third_pairs.dropna(), all_pairs.dropna()
            else:
                nn_prop, nn2_prop, nn3_prop = cat[str(props_dict[prop][0])], cat[str(props_dict[prop][1])], cat[str(props_dict[prop][2])]
                first_pairs = pd.DataFrame({'cloud_index':cloud-1, 'nn_index':nn_index, 'nn_'+str(prop):nn_prop})
                second_pairs = pd.DataFrame({'cloud_index':cloud-1, 'nn2_index':nn2_index, 'nn2_'+str(prop):nn2_prop})
                third_pairs = pd.DataFrame({'cloud_index':cloud-1, 'nn3_index':nn3_index, 'nn3_'+str(prop):nn3_prop})
                all_pairs = pd.DataFrame({'cloud_index':cloud-1, 'nn_index':nn_index, 'nn2_index':nn2_index, 'nn3_index':nn3_index,
                'nn_'+str(prop):nn_prop, 'nn2_'+str(prop):nn2_prop, 'nn3_'+str(prop):nn3_prop})
                return first_pairs.dropna(), second_pairs.dropna(), third_pairs.dropna(), all_pairs.dropna()


    if noise=='matched':
        fp = loc+str(source)+'/matched/'+str(source)+'_12m+7m+tp_co21_'+str(res)+'pc_props.fits.bz2'
        if os.path.isfile(fp):
            cat = pd.read_csv(fp)
            tab = Table.read(loc+str(source)+'/'+str(source)+'_12m+7m+tp_co21_'+str(res)+'pc_props.fits.bz2')
            cloud = cat('cloudnum')
            nn_index, nn2_index, nn3_index = cat('nn_index'), cat('nn2_index'), cat('nn3_index')
            if prop == 'MLUM_MSUN' or prop=='SIGV_KMS' or prop=='RAD_PC':
                cloud_prop = tab[prop]
                nn_prop, nn2_prop, nn3_prop = cat[str(props_dict[prop][0])], cat[str(props_dict[prop][1])], cat[str(props_dict[prop][2])]
                first_pairs = pd.DataFrame({'cloud_index':cloud-1, 'nn_index':nn_index, 'cloud_prop':cloud_prop, 'nn_'+str(prop):nn_prop})
                second_pairs = pd.DataFrame({'cloud_index':cloud-1, 'nn2_index':nn2_index, 'cloud_prop':cloud_prop, 'nn2_'+str(prop):nn2_prop})
                third_pairs = pd.DataFrame({'cloud_index':cloud-1, 'nn3_index':nn3_index, 'cloud_prop':cloud_prop, 'nn3_'+str(prop):nn3_prop})
                all_pairs = pd.DataFrame({'cloud_index':cloud-1, 'nn_index':nn_index, 'nn2_index':nn2_index, 'nn3_index':nn3_index,
                'cloud_prop':cloud_prop, 'nn_'+str(prop):nn_prop, 'nn2_'+str(prop):nn2_prop, 'nn3_'+str(prop):nn3_prop})
                return first_pairs.dropna(), second_pairs.dropna(), third_pairs.dropna(), all_pairs.dropna()
            else:
                nn_prop, nn2_prop, nn3_prop = cat[str(props_dict[prop][0])], cat[str(props_dict[prop][1])], cat[str(props_dict[prop][2])]
                first_pairs = pd.DataFrame({'cloud_index':cloud-1, 'nn_index':nn_index, 'nn_'+str(prop):nn_prop})
                second_pairs = pd.DataFrame({'cloud_index':cloud-1, 'nn2_index':nn2_index, 'nn2_'+str(prop):nn2_prop})
                third_pairs = pd.DataFrame({'cloud_index':cloud-1, 'nn3_index':nn3_index, 'nn3_'+str(prop):nn3_prop})
                all_pairs = pd.DataFrame({'cloud_index':cloud-1, 'nn_index':nn_index, 'nn2_index':nn2_index, 'nn3_index':nn3_index,
                'nn_'+str(prop):nn_prop, 'nn2_'+str(prop):nn2_prop, 'nn3_'+str(prop):nn3_prop})
                return first_pairs.dropna(), second_pairs.dropna(), third_pairs.dropna(), all_pairs.dropna()



def rad_corr(source, res, noise, prop, rad_bins, show_plot):
    global radregion
    pairs = pull_pairs(source, res, noise, prop)
    fp = '/Users/josh/projects/intro/'+str(source)+'/'+str(source)+'_12m+7m+tp_co21_'+str(res)+'pc_props.fits.bz2'
    if os.path.isfile(fp):
        tab = Table.read(fp)
        ref = np.array(pairs[3]['cloud_index'])
        pairs_in_bin = np.zeros((len(pairs[3]), pairs[3].shape[1]))
        pair_list = pairs[3].values.tolist()
        radregion = []
        for i in range(len(rad_bins)):
            for j in range(len(pairs[3])):
                if tab['RGAL_KPC'][ref[j]] > rad_bins[i-1] and tab['RGAL_KPC'][ref[j]] < rad_bins[i]:
                    #^ this line checks the index of each cloud to see if it lies within a radial bin.
                    # If it DOES, then the properties of the first, second and third nearest neighbors are retrieved
                    pairs_in_bin[j] = pair_list[j]
            radregion.append(pairs_in_bin[~(pairs_in_bin==0).all(1)])
        radregion = radregion[1:len(radregion)]
        if show_plot == True:
            fig_size = (5,10)
            fig1, ax1 = plt.subplots(len(radregion), 1, figsize=fig_size) #CLOUD PROP VS NN PROP
            fig2, ax2 = plt.subplots(len(radregion), 1, figsize=fig_size) #1st NN PROP VS 2nd NN PROP
            fig3, ax3 = plt.subplots(len(radregion), 1, figsize=fig_size) #2nd NN PROP VS 3rd NN PROP
            fig4, ax4 = plt.subplots(len(radregion), 1, figsize=fig_size) #CLOUD PROP VS 3rd NN PROP
            if i!=0:
                for k in range(len(radregion)):

                    #PLOT CLOUD PROP v 1st NEIGHBOR PROP
                    ax1[k].scatter(radregion[k][:,4], radregion[k][:,5], alpha=0.3, label='N = '+str(int(len(radregion[k][:,4])))+', '+str('%.2f' % rad_bins[k])+' < RGAL < '+str('%.2f' % rad_bins[k+1])+r'kpc, $\rho$ = '+str('%.2f' % scipy.stats.spearmanr(radregion[k][:,4], radregion[k][:,5], nan_policy='omit')[0]))
                    ax1[k].legend()
                    ax1[k].get_xaxis().set_visible(False)

                    #PLOT 1ST NN PROP v 2ND NEIGHBOR PROP
                    ax2[k].scatter(radregion[k][:,5], radregion[k][:,6], alpha=0.3, label='N = '+str(int(len(radregion[k][:,5])))+', '+str('%.2f' % rad_bins[k])+' < RGAL < '+str('%.2f' % rad_bins[k+1])+r'kpc, $\rho$ = '+str('%.2f' % scipy.stats.spearmanr(radregion[k][:,5], radregion[k][:,6], nan_policy='omit')[0]))
                    ax2[k].legend()
                    ax2[k].get_xaxis().set_visible(False)

                    #PLOT 2ND NN PROP v 3RD NEIGHBOR PROP
                    ax3[k].scatter(radregion[k][:,6], radregion[k][:,7], alpha=0.3, label='N = '+str(int(len(radregion[k][:,6])))+', '+str('%.2f' % rad_bins[k])+' < RGAL < '+str('%.2f' % rad_bins[k+1])+r'kpc, $\rho$ = '+str('%.2f' % scipy.stats.spearmanr(radregion[k][:,6], radregion[k][:,7], nan_policy='omit')[0]))
                    ax3[k].legend()
                    ax3[k].get_xaxis().set_visible(False)

                    #PLOT CLOUD PROP v 3RD NEIGHBOR PROP
                    ax4[k].scatter(radregion[k][:,4], radregion[k][:,7], alpha=0.3, label='N = '+str(int(len(radregion[k][:,4])))+', '+str('%.2f' % rad_bins[k])+' < RGAL < '+str('%.2f' % rad_bins[k+1])+r'kpc, $\rho$ = '+str('%.2f' % scipy.stats.spearmanr(radregion[k][:,4], radregion[k][:,7], nan_policy='omit')[0]))
                    ax4[k].legend()
                    ax4[k].get_xaxis().set_visible(False)

        ax1[0].set_title('Cloud vs Nearest Neighbor')
        ax1[k].set_xlabel('Cloud '+prop)
        ax1[k].get_xaxis().set_visible(True)
        fig1.text(0.01,0.5, 'Nearest Neighbor '+prop, va='center', rotation='vertical')
        fig1.suptitle(prop+' in '+str(source).upper()+' '+str(res)+'pc Resolution', fontsize=16)

        ax2[0].set_title('Nearest Neighbor vs. 2nd Nearest Neighbor')
        ax2[k].set_xlabel('Nearest Neighbor '+prop)
        ax2[k].get_xaxis().set_visible(True)
        fig2.text(0.01,0.5, '2nd Nearest Neighbor '+prop, va='center', rotation='vertical')
        fig2.suptitle(prop+' in '+str(source).upper()+' '+str(res)+'pc Resolution', fontsize=16)

        ax3[0].set_title('2nd vs 3rd Nearest Neighbor')
        ax3[k].set_xlabel('2nd Nearest Neighbor '+prop)
        ax3[k].get_xaxis().set_visible(True)
        fig3.text(0.01,0.5, '3rd Nearest Neighbor '+prop, va='center', rotation='vertical')
        fig3.suptitle(prop+' in '+str(source).upper()+' '+str(res)+'pc Resolution', fontsize=16)

        ax4[0].set_title('Cloud vs 3rd Nearest Neighbor')
        ax4[k].set_xlabel('Cloud '+prop)
        ax4[k].get_xaxis().set_visible(True)
        fig4.text(0.01,0.5, '3rd Nearest Neighbor '+prop, va='center', rotation='vertical')
        fig4.suptitle(prop+' in '+str(source).upper()+' '+str(res)+'pc Resolution', fontsize=16)

        plt.show()
        plt.close()

        return radregion

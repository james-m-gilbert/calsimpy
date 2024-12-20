# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 13:53:16 2019

@author: jmgilbert
"""
import AuxFunctions as af
import pandas as pnd
import numpy as np
import datetime as dt
#%%

import matplotlib.pyplot as plt

import seaborn as sns
sns.set_style('ticks')


import matplotlib.ticker as tkr
import matplotlib.dates as mdates
import matplotlib.patches as mpatches
import matplotlib.lines as mlines


gw_color = ['#c8cace']  #['#919cad']
sw_colors = sns.color_palette("Blues",6)
sw_colors = sns.light_palette((150/256., 59/256., 0/256.),5)
cm = sw_colors+ gw_color

def func(x, pos):  # formatter function takes tick label and tick position
# taken from http://stackoverflow.com/questions/25973581/how-do-i-format-axis-number-format-to-thousands-with-a-comma-in-matplotlib
    s = '%d' % x
    groups = []
    while s and s[-1].isdigit():
        groups.append(s[-3:])
        s = s[:-3]
    return s + ','.join(reversed(groups))

y_format = tkr.FuncFormatter(func)  # make formatter
years = mdates.YearLocator()
decades = mdates.YearLocator(10)

years_fmt = mdates.DateFormatter('%Y')

#%%  plotting a summary of porosity, field capacity, and root depth
#    -> from IDC - Parameter inputs

def make_rootzone(rootdep, fc, ef, maxdep=8, maxef=0.5):
    fig, ax = plt.subplots(1,1,figsize=(4,6))
    
    #maxdep = 8
    #maxef = 0.5
    #rootdep =  6
    #fc = 0.19
    #ef = 0.33
    
    bg = plt.Rectangle((0,0), width=maxef, height=maxdep, fill=False,ec='0.8', hatch='//', zorder=0)
    rz = plt.Rectangle((0,0), width=ef, height=rootdep, fc='lightblue', alpha=1)
    fcz = plt.Rectangle((0,0), width=fc, height=rootdep, fc='darkblue', alpha=0.8)
    #rzline = plt.Line2D(xdata=(-0.2, ef), ydata=(rootdep,rootdep),linewidth=1)
    ax.add_patch(bg)
    ax.add_patch(rz)
    ax.add_patch(fcz)
    #ax.add_line(rzline)
    ax.set_ylim((maxdep,0))
    ax.set_xlim((0, maxef))
    ax.xaxis.tick_top()
    #ax.set_ylabel('Depth from ground surface, ft')
    
    transx = ax.get_xaxis_transform()
    ax.annotate('Root\ndepth', xy=(0, rootdep), xycoords='data', xytext=(-65, -5), textcoords='offset points', #xy=(ef, rootdep),xycoords=transx, xytext=(-0.5,rootdep),
                fontstyle='italic', fontsize=10) #, arrowprops=dict(arrowstyle="-"))
    ax.plot([-0.09, -0.05], [rootdep, rootdep],lw=0.8,ls='-', color="k", clip_on=False) #, transform=transx, clip_on=False)
    
    
    ax.annotate('Porosity', xy=(ef, 0), xycoords='data', xytext=(-8, 40), textcoords='offset points', #xy=(ef, rootdep),xycoords=transx, xytext=(-0.5,rootdep),
                fontstyle='italic', fontsize=10) #, arrowprops=dict(arrowstyle="-"))
    ax.plot([ef, ef], [rootdep, -0.5],lw=1.5,ls='-', color="gray", clip_on=False)
    
    
    ax.annotate('Field\ncapacity', xy=(fc, 0), xycoords='data', xytext=(-8, 30), textcoords='offset points', #xy=(ef, rootdep),xycoords=transx, xytext=(-0.5,rootdep),
                fontstyle='italic', fontsize=10) #, arrowprops=dict(arrowstyle="-"))
    ax.plot([fc, fc], [rootdep, -0.5],lw=1.5,ls='-', color="gray", clip_on=False)
    
    ax.set_ylabel('Depth, ft')
    ax.set_xlabel('Volume fraction')
    plt.subplots_adjust(left=0.3, top=0.85)
    
    return([fig,ax])

def patternizer(tsdf, freq='M', **kwargs):
    
    if type(tsdf) ==pnd.Series:
        tt1 = pnd.DataFrame(data=tsdf.copy(), columns=['value'])
    else:
        tt1 = tsdf.copy()      

    
    if freq=='M':    
        tt1['Month'] = tt1.index.month
        tt1['WYmonthorder'] = tt1.index.map(lambda x: af.wymo(x))
        tt1['Year'] = tt1.index.year
        tt1['WY'] = tt1.index.map(lambda x: af.addWY(x))
        
    if 'reference' in kwargs:
        if kwargs['reference']=='WY':    
            ttpv = pnd.pivot(tt1, index='WYmonthorder', columns='WY')[tt1.columns[0]] #['value']
            ttpv['Month'] = [10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    
    return(ttpv)
    

def patternPlot(patDF, indexCol='', xLabelCol='Month',xLabelFmt='%b'):
    
    df = patDF.copy()
    idxname = df.index.name
    df.reset_index(inplace=True)
    
    fig, ax = plt.subplots(1,1,figsize=(6,4))
    
    if indexCol =='index':
        idx = df.pop(idxname)
        
    elif indexCol not in df:
        print("The prescribed indexCol variable %s is not in the dataframe that was passed" %indexCol)
        return
    else: 
        idx = df.pop(indexCol)
        
    idxlab = []
    for i in df[xLabelCol]:
        tmpdt = dt.datetime(2020, i, 1, 0, 0)
        idxlab.append(tmpdt.strftime(xLabelFmt))
        
    
    for c in df.columns:
        if type(c)==str and c.upper() not in ['MONTH','YEAR', 'DAY']:
            ax.plot(idx, df.loc[:,c], color='0.8',lw=0.8)
        elif type(c)!=str:
            ax.plot(idx, df.loc[:,c], color='0.8',lw=0.8)
        else:
            pass
    ax.set_xticks(idx)
    ax.set_xticklabels(idxlab)    
    return([fig,ax])
    


#%% Plot the total applied water and sources (GP, net deliveries)  
def plot_AppliedWater_Annual(duID, csDF, cshDF, contract_type, year_agg_method=2):    
    
    '''
        assumes csDF is a dataframe of calsim delivery results with columns
        indicating the applied demand, net& gross deliveries, & gw pumping, 
        riparian use, reuse, etc; all in volumes TAF
        
        cshDF is a datafream of calsimhydro budget data, all in volumes of TAF
        
        
        # year_agg_method: 1 = Water year (Oct-Sep); 2 = Contract year (Mar-Feb)
        
    '''
   
    plotcols =1 

    if duID=='71_PA4':  #TODO: this is a hack to deal with du's that might have multiple delivery-contract type combinations
        csDF['clmVar'] = csDF['clmVar_a'] + csDF['clmVar_b']
        csDF['clmAllocVar'] = csDF['clmVar_a'] * csDF['allocAg'] + csDF['clmVar_b'] * csDF['allocEX']
#        if 'allocVar' in csDF.columns: # calculate the actual allocation contract limit
#            csDF['clmAllocVar'] = csDF['clmVar'] * csDF['allocVar']
#        else:
#            csDF['clmAllocVar'] = np.nan*len(csDF['clmVar'])
    else:
        if contract_type.upper()=='AG':
            allocCol='allocAg'
        elif contract_type.upper()=='EX':
            allocCol='allocEX'
        elif contract_type.upper()=='REF':
            allocCol='allocRef'
        elif contract_type.upper()=='MI':
            allocCol='allocMI'
        else:
            allocCol = 'allocAg'
        csDF['clmAllocVar'] = csDF['clmVar'] * csDF[allocCol]
    
    if year_agg_method == 1:
        csh_ann = cshDF.resample('A-Sep').apply(np.sum)       
        cs3_ann = csDF.resample('A-Sep').apply(np.sum)  
    elif year_agg_method == 2:
        csh_ann = cshDF.resample('A-Feb').apply(np.sum)
        cs3_ann = csDF.resample('A-Feb').apply(np.sum)      
    else:  # if not WY or Contract year, assume calendar year (Jan-Dec)
        csh_ann = cshDF.resample('A-Dec').apply(np.sum)
        cs3_ann = csDF.loc[:,:].resample('A-Dec').apply(np.sum) 
    
    with sns.plotting_context('paper', font_scale=1.5):
        fig, axes = plt.subplots(1,plotcols,figsize=(11,6))    #<--line to use if doing regular sublots
        #fig = plt.figure(figsize=(12,13))
        #grid = plt.GridSpec(3,3, wspace=0.4, hspace=0.3)   
        #area=df_idc.loc[1,:][du+'_'+landuse+'_' +'AREA'].iloc[0]
        
        #tmp_idc_ann = csh_ann['TOTAL_APP']
        tmp_idc_ann = cs3_ann['awVar']
        
        if 'NA' in duID or 'PA' in duID or 'PR' in duID or 'XA' in duID:  #ag lands
            ru_idx = [i for i, s in enumerate(cs3_ann.columns) if 'RU' in s.upper()][0]
            rp_idx = [i for i, s in enumerate(cs3_ann.columns) if 'RP' in s.upper()][0]
            tmp_idc_ann = tmp_idc_ann - cs3_ann.iloc[:,ru_idx]+cs3_ann.iloc[:,rp_idx]
        elif 'NU' in duID or 'PU' in duID:  #urban
            tmp_idc_ann = tmp_idc_ann
        else:
            tmp_idc_ann = tmp_idc_ann
            
            
        #axes = tmp_idc_TAF.resample('A-Sep').apply(np.sum).plot(label='Total Applied Water', lw=3, c='w',drawstyle="steps-post" )

            
        axes.plot(tmp_idc_ann.index, tmp_idc_ann,
                  label='Total Applied Water Demand (AWR+AWO)\nPlus Riparian/Misc ET\nMinus Reuse', 
                  lw=2, c='orange',drawstyle="steps-post" )
        
        numHndles=2
        if 'PA' in duID or 'PU' in duID or 'XA' in duID or 'PR' in duID:
            numHndles = 3
            axes.plot(cs3_ann.index, cs3_ann.clmVar,
                      label='Contract Limit', 
                      lw=2, c='k',drawstyle="steps-post" )
            axes.plot(cs3_ann.index, cs3_ann.clmAllocVar,
                      label='Allocation limit',
                      lw=0.7, c='k', ls='-', drawstyle="steps-post")
    
        #cs3_df_ann.plot(kind='area', ax=axes, alpha=0.6, drawstyle="steps-post")
        dnCols = [s for i, s in enumerate(cs3_ann.columns) if 'DN' in s.upper()]
        if len(dnCols)==1:
            dnColNames = ['Net Delivery']
            sw_colors = sns.light_palette((20/256., 59/256., 150/256.),5)
            sw_colors.reverse()
            cm = [sw_colors[1]]+ gw_color
        else:
            dnColNames = ['NetDel - %s' %x[0:10] for x in dnCols ]
            sw_colors = sns.light_palette((20/256., 59/256., 150/256.),min(5, len(dnCols)+1))
            sw_colors.reverse()
            cm = sw_colors[1:] + gw_color
        
        gpCols = [s for i, s in enumerate(cs3_ann.columns) if 'GP' in s.upper()]
        
        if len(gpCols)==1:
            gpColNames = ['GW Pumping']
        else:
            gpColNames = ['GW Pumping - %s' %x[0:10] for x in gpCols ]
    
        iterCols =  dnCols + gpCols
        iterColNames = dnColNames + gpColNames
        prev = 0
        for n,c in enumerate(zip(iterCols, iterColNames)): #cs3_ann.columns[:-2]):
            axes.fill_between(cs3_ann.index, cs3_ann[c[0]]+prev, prev, 
                              lw=0.,step='post', color=cm[n],alpha=0.99,
                              label =c[1])    #'_'.join(c.split('_')[0:-2]) )
            prev = cs3_ann[c[0]]+prev
    
        handles, labels = axes.get_legend_handles_labels()
        handles2 = []
        labels2 = []
    
        print(labels)
        for l in labels:
            l2 = l 
            labels2.append(l2)
    
        for ih in range(numHndles):
            handles2.append(handles[ih])
            #handles2.append(handles[1])
        for h in handles[numHndles:]:
            c1 = h.get_facecolor()
            #h2 = mlines.Line2D([],[], color=c1, linewidth=4)
            h2 = mpatches.Patch( facecolor=c1[0])
            handles2.append(h2)
    
        plt.subplots_adjust(bottom=0.21, top=0.93, left=0.1, right=0.95)
    
        axes.set_xlim(dt.date(1920,1,1), dt.date(2020, 1,1))
        axes.xaxis.set_major_locator(decades)
        axes.xaxis.set_major_formatter(years_fmt)
        

        axes.xaxis.set_minor_locator(years)
        
        if year_agg_method==1:
            axes.set_ylabel('Water Year Volume (Oct-Sep), TAF')
        elif year_agg_method==2:
            axes.set_ylabel('Contract Year Volume (Mar-Feb), TAF')
        else:
            axes.set_ylabel('Calendar Year Volume (Jan-Dec), TAF')
        
        if axes.get_ylim()[1]>100:
            axes.yaxis.set_major_locator(tkr.MultipleLocator(100))
            axes.yaxis.set_minor_locator(tkr.MultipleLocator(20))
            axes.yaxis.set_major_formatter(y_format)
            
        elif axes.get_ylim()[1]<10:
            axes.yaxis.set_major_locator(tkr.MultipleLocator(1))
            axes.yaxis.set_minor_locator(tkr.MultipleLocator(0.2))
            
        elif axes.get_ylim()[1]<50:
            axes.yaxis.set_major_locator(tkr.MultipleLocator(5))
            axes.yaxis.set_minor_locator(tkr.MultipleLocator(1))
            
        elif axes.get_ylim()[1]<100:
            axes.yaxis.set_major_locator(tkr.MultipleLocator(20))
            axes.yaxis.set_minor_locator(tkr.MultipleLocator(5))
            axes.yaxis.set_major_formatter(y_format)
        else:
            axes.yaxis.set_major_locator(tkr.MultipleLocator(100))
            axes.yaxis.set_minor_locator(tkr.MultipleLocator(20))
            axes.yaxis.set_major_formatter(y_format)
    
        
        axes.set_title('Demand Unit: ' + duID)
        leg1 = fig.legend(handles=handles2, labels=labels2, loc='lower center',
                          bbox_to_anchor=(0.5,-0.01), ncol=4, fontsize=11,
                          frameon=False)
        #leg2 = fig.legend(handles=handles_neg, labels=labels_neg, loc='lower center',bbox_to_anchor=(0.75,-0.01), ncol=2, fontsize=12 )
        #plt.tight_layout()

        sns.despine()
        
        return([fig,axes])

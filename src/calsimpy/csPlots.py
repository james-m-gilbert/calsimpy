# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 15:04:17 2020

@author: jgilbert
"""

from collections import OrderedDict as Odict
import cs3
import AuxFunctions as af
import os, sys
import datetime as dt
import numpy as np

from pathlib import Path

import pandas as pnd
idx = pnd.IndexSlice

import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('ticks')

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoLocator,
                               AutoMinorLocator, PercentFormatter)

import matplotlib.dates as mdates

from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

from copy import deepcopy as dcp

years = mdates.YearLocator()   # every year
months = mdates.MonthLocator()  # every month
years_fmt = mdates.DateFormatter('%Y')

decades = {'1920-1940': ['1920-10-31','1939-09-30'],
           '1940-1960': ['1939-10-31','1959-09-30'],
           '1960-1980': ['1959-10-31', '1979-09-30'],
           '1980-2000': ['1979-10-31', '1999-09-30'],
           '2000-2020': ['1999-10-31', '2020-09-30']}

month_dict= {1:'Jan',2:'Feb',3:'Mar',4:'Apr',5:'May',6:'Jun',7:'Jul',
             8:'Aug',9:'Sep',10:'Oct',11:'Nov',12:'Dec'}


def af_to_mcm(afv):
    mcmv = afv/1e6*43560.*(12**3)*(2.54**3)/(100**3)
    return(mcmv)
def mcm_to_af(mcmv):
    mafv = mcmv*(100**3)/(2.54**3)/(12**3)/43560
    return(mafv*1e6)

def taf_to_mcm(tafv):
    mcmv = (tafv*1000.)/1e6*43560.*(12**3)*(2.54**3)/(100**3)
    return(mcmv)
def mcm_to_taf(mcmv):
    mafv = mcmv*(100**3)/(2.54**3)/(12**3)/43560
    return(mafv*1e3)

class cs_analysis:
    
    def __init__(self,csStudies):
        
        self.Studies = Odict()
        if type(csStudies)==list:
            for i, c in enumerate(csStudies):
                self.Studies[i] = c
        else:
            idx = max(self.Studies.keys(),0)
            self.Studies[idx+1] = csStudies
            
#        if 'study_dict' in kwargs:
#            self.Studies= kwargs['study_dict']
#        else:
#            self.Studies = {}
        
        self.Analytes= Odict()
        self.Staged = Odict()
        
    def include_study(self, csStudy):  #csObj, color, short_name, desc='', baseline=False):
        idx = max(self.Studies.keys(),0)
        self.Studies[idx+1] = csStudy
        
        
    def getAllSV(self):
        if len(self.Studies)>0:
            for s in self.Studies:
                self.Studies[s].getSV()
                

    def getDV(self, filter=''):
        '''
            convenience function to get DV time series
            `filter` can be a single string or list, as long as each string
                     is a valid contiguous set of characters that could be found
                     in a DSS record path (e.g. '/C5/','/FLOW-INFLOW')
                     A blank string provided to filter (default) will retrieve
                     all DV ts data
        '''
        if len(self.Studies)>0:
            for s in self.Studies:
                self.Studies[s].getDV(filter)
                
    def getSV(self, filter=''):
        '''
            convenience function to get SV time series
            `filter` can be a single string or list, as long as each string
                     is a valid contiguous set of characters that could be found
                     in a DSS record path (e.g. '/C5/','/FLOW-INFLOW')
                     A blank string provided to filter (default) will retrieve
                     all DV ts data
        '''
        if len(self.Studies)>0:
            for s in self.Studies:
                self.Studies[s].getSV(filter)

    def getAnalyte(self,analysisVar, SVDV='SV'):
        '''
            assuming analysisVar is a single B-part variable name for now
            TODO: add multivariable/derived var functionality
        '''
        if len(self.Studies)>0:
            tmp = Odict()
            for s in self.Studies:
                if SVDV == 'SV':
                    tmp[s] = self.Studies[s].CalSimObj.SVdata.SVtsDF.loc[:,idx[:,analysisVar]]
                else:
                    tmp[s] = self.Studies[s].CalSimObj.DVdata.DVtsDF.loc[:,idx[:,analysisVar]]
            self.Analytes[analysisVar] = tmp
                
    def stageVar(self, var, SVDV,**kwargs):
        '''
            an alternative to `getAnalyte`, meant to be called 
            by the plotting functions
            
            still assuming var is a single B-part variable name for now
            TODO: add multivariable/derived var functionality
            
        '''
        
        if 'varname' in kwargs:
            varname = kwargs['varname']
        else:
            varname = var
            
        if len(self.Studies)>0:
            tmp = Odict()
            tmp1 = Odict()
            for s in self.Studies:
                
                if SVDV == 'SV':
                    if type(var)==list:
                        for v in var:
                            tmp1[v] = self.Studies[s].CalSimObj.SVdata.SVtsDF.loc[:,idx[:,v]]
                        tmp2 = pnd.concat(tmp1, axis=1)
                        # assuming addition for now - would need a way to parse
                        # the operators
                        tmp[s] = tmp2.sum(axis=1)
                    else:
                        tmp[s] = self.Studies[s].CalSimObj.SVdata.SVtsDF.loc[:,idx[:,var]]
                else:
                    if type(var)==list:
                        for v in var:
                            tmp1[v] = self.Studies[s].CalSimObj.DVdata.DVtsDF.loc[:,idx[:,v]]
                        tmp2 = pnd.concat(tmp1, axis=1)
                        bcol = tmp2.columns[0]
                        # assuming addition for now - would need a way to parse
                        # the operators
                        tmp3 = pnd.DataFrame(tmp2.sum(axis=1))
                        tmp3.columns = [['CALSIM'],[varname],[bcol[3]],[bcol[4]],[bcol[5]],[bcol[6]],[bcol[7]]]
                        cidx = tmp3.columns.set_names(['A','B','C','E','F','Type','Units'])
                        tmp3.columns = cidx
                        tmp[s] = tmp3
                    else:
                        tmp[s] = self.Studies[s].CalSimObj.DVdata.DVtsDF.loc[:,idx[:,var]]
            
            self.Staged[varname] = tmp
            
            
    def unstageVar(self, var):
        delVar = self.Staged.pop(var, None)
        print("Un-staged variable: %s" %delVar.keys())
        #self.Staged = self.Staged
    
    def stageList(self):
        for s in self.Staged:
            print(s)
    
    
    def plot_annual_exceed_shaded_pctile(self, title, annualize_on='A-Sep', how='auto-sum',
                           month_filter=[],annual_filter=[],
                           annual_filter_type="WY",
                           annual_filter_label="", **kwargs):

        '''
            month_filter: will plot exceedance based on just months meeting filter criteria (e.g. month_filter=[9] will include only September)
            annual filter: will plot exceedance based on just years meeting filter criteria
            
            
            kwargs:
                'reverse_x': reverses x-axis (exceedance probabilities) so that 100% is at the left
        '''
        fig_dict = {}
        
        if 'exclude' in kwargs:
            exclude_scen = kwargs['exclude']
        else:
            exclude_scen = []
            
        if type(title)!=list:
            title = [title]
            
        var_cntr = 0
        for var, dictDF in self.Staged.items():   # iterate through the staged variables, each one having a set of dataframes
        
            numpctiles = 11
            from matplotlib import cm
            colormap = cm.Blues
            
            pctiles = np.linspace(0, 100, numpctiles)
            
            #master_df = pnd.DataFrame(index=dictDF[list(dictDF.keys())[0]].index)
            master_list = []
            nscens = len(dictDF) # number of scenarios
            cntr = 0
            for scen, df in dictDF.items():   # iterate through the scenarios and dataframes
            
                if scen in exclude_scen:
                    continue
                
                if self.Studies[scen].Baseline:
                    bl_scen = scen
                    isBaseline=True
                    thiszorder = 9999
                else:
                    isBaseline=False
                    thiszorder = 0
                
                origUnits = df.columns.get_level_values('Units')[0]
                
                dftmp = df.copy()
                dftmp['WY'] = dftmp.index.map(af.addWY)
                
                dftmp = dftmp.loc['1921-10-31':'2021-09-30',:]
                # FIRST: filter out jsut the years or months desired
                
                # check for annual filters - this filters out specific years
                if annual_filter == []:
                    annfilt = dftmp
                else:
                    if annual_filter_type.upper() in ["WY","WATERYEAR","WATER YEAR"]:
                        annfilt = dftmp[dftmp.WY.isin(annual_filter)]
                    else:
                        annfilt = dftmp[dftmp.index.year.isin(annual_filter)]  #??? This only filters on calendar year - 
                                                                     # TODO: will need to add a way to do water year or delivery year
                                                                     # type filters
                
                # check for month filters - this selects for specific months
                if month_filter == []:
                    monfilt = annfilt
                else:
                    monfilt = annfilt[annfilt.index.month.isin(month_filter)]
                    
                # drop the WY column to avoid the hassle of carrying that around
                #print(monfilt.head())
                filt = monfilt.drop(columns=['WY'])
                
                # SECOND: now do the aggregation/annualization
                if annualize_on != None:
                    agg_df= annualize(filt, on=annualize_on, how=how, 
                                      colindex=[0], )
                    finUnits = agg_df.columns.get_level_values('Units')[0]
                else:
                    agg_df = filt
                    finUnits = origUnits
                            
                excd_df = single_exceed(agg_df, agg_df.columns[0])
                
                if isBaseline:
                    baseline_df = excd_df
                else:
                    if cntr==0:
                        master_df = pnd.DataFrame(index=excd_df.index)
                        master_df[scen] = excd_df
    
                    else:
                        master_df[scen] =  excd_df
                    #master_list.append(excd_df.iloc[0:].values)
                    cntr+=1
                
            #master_df = pnd.DataFrame(data=master_list, index=excd_df.index)

            sdist = np.zeros((len(master_df), numpctiles))
            for i in range(numpctiles):
                for t,d in enumerate(master_df.index):
                    sdist[t, i] = np.percentile(master_df.loc[d,:], pctiles[i])
            maxt = sdist[:,-1]
            mint = sdist[:, 0]
            half = int((numpctiles-1)/2)

            with sns.plotting_context('notebook', font_scale=1.5):
                fig, ax = plt.subplots(1,1, figsize=(11,8))
                ax.plot(master_df.index, sdist[:, half], c='k')
                ax.plot(master_df.index, np.mean(master_df, axis=1), ls='--',lw=2, c='k')
                ax.plot(baseline_df.index, baseline_df.values,c='darkred', lw=2)
                for i in range(half):
                    ax.fill_between(master_df.index, sdist[:,i], sdist[:,-(i+1)], 
                                    color=colormap(i/half))
                ax.plot(master_df.index, mint, lw=0.5, c='k')
                ax.plot(master_df.index, maxt, lw=0.5, c='k')
        
                ax.set_xlabel('Exceedance probability')
                
                if 'reverse_x' in kwargs:
                    if kwargs['reverse_x']:
                        ax.invert_xaxis()
                        
                ax.xaxis.set_major_formatter(PercentFormatter(1.0))
                ax.xaxis.set_major_locator(MultipleLocator(0.1))
                ax.xaxis.set_minor_locator(AutoMinorLocator())
                ax.yaxis.set_minor_locator(AutoMinorLocator())
                
                if finUnits.upper() in ['TAF','AF','ACRE-FEET']:
                    ax.set_ylabel('Annual Volume (%s)' %finUnits)
                else:
                    ax.set_ylabel('Annual Average Flow (%s)' %finUnits)   #TODO: this is just a placeholder - need to add ability to deal with cfs, EC, other units and aggregations
                
                if month_filter != []:
                    monthTitle = "months: "
                    for m in month_filter:
                        if len(month_filter)==1:
                            monthTitle = f"{dt.date(2020,m,1):%B}"
                        else:
                            monthTitle += "%s, " %month_dict[m]
                else:
                    monthTitle = ""
                
                if annual_filter != [] and annual_filter_label !="":
                    annTitle = annual_filter_label
                elif annual_filter != []:
                    annTitle = "Select years"
                else:
                    annTitle = ""
                    
                if annTitle=="":
                    sep=''
                else:
                    sep='-'
                full_title = title[var_cntr] + " Annual Exceedance\n%s %s %s" %(annTitle,sep, monthTitle)
                ax.set_title(full_title, fontsize=14, fontweight='bold')
            
            var_cntr+=1
                
            # if nscens > 7:
            #     h1, l1 = ax.get_legend_handles_labels()
            #     h2 = [h1[0]]+[h1[1]] # get baseline + first realization handles
            #     l2 = [l1[0]] + ['realizations']
            #     plt.legend(h2,l2, bbox_to_anchor=(0.2,-0.25, 0.6, 0.15), ncol=2,
            #                frameon=False, fontsize=11)
            # else:
            #     legcols = nscens%2 + nscens//2
            #     plt.legend(bbox_to_anchor=(0.2,-0.25, 0.6, 0.15), ncol=legcols,frameon=False, fontsize=11)
                
            plt.subplots_adjust(bottom=0.22,top=0.90) #,hspace=0.55)
            sns.despine()
            fig_dict[var] = [fig, ax]        

        return(fig_dict)
        
    def plot_annual_exceed(self,title, annualize_on='A-Sep', how='auto-sum',
                           month_filter=[],annual_filter=[],
                           annual_filter_type="WY",
                           annual_filter_label="", 
                           return_data = False,
                           **kwargs):
        '''
            month_filter: will plot exceedance based on just months meeting filter criteria (e.g. month_filter=[9] will include only September)
            annual filter: will plot exceedance based on just years meeting filter criteria
            
            
            kwargs:
                'reverse_x': reverses x-axis (exceedance probabilities) so that 100% is at the left
        '''
        fig_dict = {}
        
        if 'exclude' in kwargs:
            exclude_scen = kwargs['exclude']
        else:
            exclude_scen = []
            
        if type(title)!=list:
            title = [title]
            
        excd_data = {}
        var_cntr = 0
        for var, dictDF in self.Staged.items():   # iterate through the staged variables, each one having a set of dataframes
            
            var_data_ = {}
            with sns.plotting_context('notebook', font_scale=1.2):
                fig, ax = plt.subplots(1,1, figsize=(10,8))
    
                nscens = len(dictDF)
                for scen, df in dictDF.items():   # iterate through the scenarios and dataframes
                    scen_name = self.Studies[scen].ShortName
                    
                    if (scen in exclude_scen) or (scen_name in exclude_scen):
                        continue
                    
                    if self.Studies[scen].Baseline:
                        isBaseline=True
                        lw=2
                        thiszorder = 9999
                    else:
                        isBaseline=False
                        if 'lw' in kwargs:
                            lw = kwargs['lw']
                        else:
                            lw=1.2
                        thiszorder = 0
                    
                    ls = self.Studies[scen].LineStyle
                    
                    origUnits = df.columns.get_level_values('Units')[0]
                    
                    dftmp = df.copy()
                    dftmp['WY'] = dftmp.index.map(af.addWY)
                    
                    dftmp = dftmp.loc['1921-10-31':'2021-09-30',:]
                    # FIRST: filter out jsut the years or months desired
                    
                    # check for annual filters - this filters out specific years
                    if annual_filter == []:
                        annfilt = dftmp
                    else:
                        if annual_filter_type.upper() in ["WY","WATERYEAR","WATER YEAR"]:
                            annfilt = dftmp[dftmp.WY.isin(annual_filter)]
                        else:
                            annfilt = dftmp[dftmp.index.year.isin(annual_filter)]  #??? This only filters on calendar year - 
                                                                         # TODO: will need to add a way to do water year or delivery year
                                                                         # type filters
                    
                    # check for month filters - this selects for specific months
                    if month_filter == []:
                        monfilt = annfilt
                    else:
                        monfilt = annfilt[annfilt.index.month.isin(month_filter)]
                        
                    # drop the WY column to avoid the hassle of carrying that around
                    #print(monfilt.head())
                    filt = monfilt.drop(columns=['WY'])
                    
                    # SECOND: now do the aggregation/annualization
                    if annualize_on != None:
                        agg_df= annualize(filt, on=annualize_on, how=how, 
                                          colindex=[0], )
                        finUnits = agg_df.columns.get_level_values('Units')[0]
                    else:
                        agg_df = filt
                        finUnits = origUnits
#                        if annual_filter == []:
#                            ann_df= annualize(dftmp, on=annualize_on, how=how, colindex=[0])
#                        else:
#                            if annual_filter_type.upper() in ["WY","WATERYEAR","WATER YEAR"]:
#                                df2 = dftmp[dftmp.WY.isin(annual_filter)]
#                            else:
#                                df2 = dftmp[dftmp.index.year.isin(annual_filter)]  #??? This only filters on calendar year - 
#                                                                         # TODO: will need to add a way to do water year or delivery year
#                                                                         # type filters
                            #ann_df = annualize(df2, on=annualize_on, how=how, colindex=[0])
                                         
#                        if month_filter==[]:
#                            ann_df = annualize(df, on=annualize_on, how=how, colindex=[0])
#                        else:
#                            df2 = df[df.index.month.isin(month_filter)]
#                            ann_df = annualize(df2, on=annualize_on, how=how, colindex=[0])
                                    
                    excd_df = single_exceed(agg_df, agg_df.columns[0])
                    
                    var_data_[scen_name] = excd_df    
                    
                    colr = self.Studies[scen].Color
                    lab = self.Studies[scen].ShortName
                    
                    if 'labels' in kwargs:
                        if kwargs['labels'].lower() in ['desc','description','descr']:
                            lab = self.Studies[scen].Description
                        else:
                            raise NotImplementedError('labels options limited to "description" for now')
                    
                    
                    if 'color' in kwargs:
                        if type(kwargs['color'])==dict:
                            if scen_name in kwargs['color']:  # by name
                                colr = kwargs['color'][scen_name]
                            elif scen in kwargs['color']: # by index
                                colr = kwargs['color'][scen]
                            if not isBaseline:
                                lw = 1
                        else:
                            if not isBaseline:
                                colr = kwargs['color']
                                
                    if 'zorder' in kwargs:
                        if type(kwargs['zorder'])==dict:
                            if scen_name in kwargs['zorder']:
                                thiszorder=kwargs['zorder'][scen_name]
                                
                    ax.plot(excd_df.index, excd_df.iloc[:,0], color= colr, 
                            label=lab, lw=lw, zorder=thiszorder, ls=ls)
                    
                    ax.set_xlabel('Exceedance probability')
                    
                    if 'reverse_x' in kwargs:
                        if kwargs['reverse_x']:
                            ax.invert_xaxis()
                            
                    ax.xaxis.set_major_formatter(PercentFormatter(1.0))
                    ax.xaxis.set_major_locator(MultipleLocator(0.1))
                    ax.xaxis.set_minor_locator(AutoMinorLocator())
                    ax.yaxis.set_minor_locator(AutoMinorLocator())
                    
                    if finUnits.upper() in ['TAF','AF','ACRE-FEET']:
                        ax.set_ylabel('Annual Volume (%s)' %finUnits)
                    else:
                        ax.set_ylabel('Annual Average Flow (%s)' %finUnits)   #TODO: this is just a placeholder - need to add ability to deal with cfs, EC, other units and aggregations
                    
                    
                    if finUnits.upper() =='TAF':
                        ax2 = ax.secondary_yaxis( -0.12, functions=(taf_to_mcm, mcm_to_taf))
                        ax2.yaxis.set_major_locator(AutoLocator()) #tkr.MultipleLocator(1))
                        ax2.yaxis.set_minor_locator(AutoMinorLocator()) #tkr.MultipleLocator(0.2))
                        ax2.set_ylabel('Million cubic meters',fontsize=14)
                        plt.subplots_adjust(left=0.18)
                    
                    if month_filter != []:
                        monthTitle = "months: "
                        for m in month_filter:
                            if len(month_filter)==1:
                                monthTitle = f"{dt.date(2020,m,1):%B}"
                            else:
                                monthTitle += "%s, " %month_dict[m]
                    else:
                        monthTitle = ""
                    
                    if annual_filter != [] and annual_filter_label !="":
                        annTitle = annual_filter_label
                    elif annual_filter != []:
                        annTitle = "Select years"
                    else:
                        annTitle = ""
                        
                    if annTitle=="":
                        sep=''
                    else:
                        sep='-'
                    if month_filter != []:
                        full_title = title[var_cntr] + " Exceedance\n%s %s %s" %(annTitle,sep, monthTitle)
                    else:
                        full_title = title[var_cntr] + " Annual Exceedance\n%s %s %s" %(annTitle,sep, monthTitle)
                    ax.set_title(full_title, fontsize=14, fontweight='bold')
            excd_data[var] = var_data_
            var_cntr+=1
                
            if nscens > 15:
                h1, l1 = ax.get_legend_handles_labels()
                #hcolors = list(set([hh.get_color() for hh in h1]))
                h2dict = {}
                for hh in h1:
                    h2dict[hh.get_color()] = hh
                h2 = [v for k,v in h2dict.items()]
                if 'group_labels' not in kwargs:
                    l2 = [l1[0]] + ['realizations']
                else:
                    l2 = kwargs['group_labels']
                #l2 = []
                #h2 = [h1[0]]+[h1[1]] # get baseline + first realization handles
                #
                plt.legend(h2,l2, bbox_to_anchor=(0.2,-0.28, 0.6, 0.15), ncol=2,
                           frameon=False, fontsize=14)
            else:
                legcols = nscens%2 + nscens//2
                #plt.legend(bbox_to_anchor=(0.3,-0.25, 0.6, 0.15), ncol=legcols,frameon=False, fontsize=11)
                axbox = ax.get_position()
                #+0.5*axbox.width
                # plt.legend(bbox_to_anchor=[axbox.x0+0.2*axbox.width, axbox.y0-0.14,0.6, 0.15], 
                #            bbox_transform=fig.transFigure,
                #            ncol=legcols,frameon=False, fontsize=11)
                fig.legend(loc='outside lower center', fontsize=15, 
                           ncol=legcols,frameon=False,)
                
            plt.subplots_adjust(bottom=0.22,top=0.90) #,hspace=0.55)
            sns.despine()
            fig_dict[var] = [fig, ax]        
#            if SCRATCH:
#                plt.savefig(os.path.join(scratch_dir, 'climate_scenarios_ShastaInflowAnn_compare.png'), dpi=300)
        if return_data:
            return([fig_dict, excd_data])
        else:
            return(fig_dict)
    


    
    def plot_pattern(self, title, ylabel, plotType='avg',
                     annual_filter=[],
                     annual_filter_type="WY",
                     annual_filter_label="",
                     pattern_type="WY",
                     indexCol='index',
                     xLabelCol='Month',
                     xLabelFmt='%b',
                     **kwargs):
        '''
            plots the monthly pattern of the time series
            TODO: filter out years by type or explicit list
        '''
    #def patternPlot(patDF, indexCol='', xLabelCol='Month',xLabelFmt='%b'):
        fig_dict = {}
        
        if 'exclude_months' in kwargs:
            exclude_months = kwargs['exclude_months']
        else:
            exclude_months = []
            
        if 'exclude' in kwargs:
            exclude_scen = kwargs['exclude']
        else:
            exclude_scen = []
            

            
        for var, dictDF in self.Staged.items():   # iterate through the staged variables, each one having a set of dataframes
            
            with sns.plotting_context('paper', font_scale=1.3):
                fig, ax = plt.subplots(1,1, figsize=(10,7))
    
                nscens = len(dictDF)
                
                scen_df = pnd.DataFrame()
                pal = {}
                labs = []
                
                nlabls = []
                nhdls = []
                
                for scen, df in dictDF.items():   # iterate through the scenarios and dataframes
                
                    if scen in exclude_scen:
                        continue
                    
                    origUnits = df.columns.get_level_values('Units')[0]
                    
                    if 'units' in kwargs:
                        newUnits = kwargs['units']
                    else:
                        newUnits = origUnits
                    
                    dftmp = df.copy()         
                    
                    if origUnits.upper() =='TAF' and newUnits.upper()=='CFS':
                        dftmp = taf_to_cfs(dftmp)
                    if origUnits.upper()=='CFS' and newUnits.upper()=='TAF':
                        dftmp = cfs_to_taf(dftmp)
                    
                    # do the 'patternizing' to get the DF set up the way we need it
                    patDF = patternizer(dftmp, reference=pattern_type)
                    data_cols = [y for y in patDF.columns if str(y).upper() not in ['CALYRMONTH','MONTH','YEAR','AVG','MAX','MIN']]
                    patDF['avg'] = patDF[data_cols].mean(axis=1, skipna=True) #patDF[select_wys].mean(axis=1)
                    patDF['max'] = patDF[data_cols].max(axis=1, skipna=True) #patDF[select_wys].max(axis=1)
                    patDF['min'] = patDF[data_cols].min(axis=1, skipna=True) #patDF[select_wys].min(axis=1)
                    
                    
                    if 'quantile' in kwargs:
                        blankrows = patDF.index[patDF[data_cols].isnull().all(axis=1)]
                        if type(kwargs['quantile'])==list:
                            for qt in kwargs['quantile']:
                                qlab = 'q%s' %(10*qt)
                                patDF[qlab] = patDF[data_cols].quantile(qt, interpolation='nearest', axis=1)
                                patDF.loc[blankrows, qlab] = np.nan
                        else:  # assume it's  single value
                            qlab = 'q%s' %(10*qt)
                            patDF[qlab] = patDF.quantile(qt, interpolation='nearest', axis=1)
                                
                    #patDF['q50'] = patDF.quantile(0.5, interpolation='nearest', axis=1)
                    #patDF['agg'] = 
                      
                    idxname = patDF.index.name
                    patDF.reset_index(inplace=True)
                                
                    if indexCol =='index':
                        idx = patDF.pop(idxname)
                        
                    elif indexCol not in patDF:
                        print("The prescribed indexCol variable %s is not in the dataframe that was passed" %indexCol)
                        return
                    else: 
                        idx = patDF.pop(indexCol)
                        
                    idxlab = []
                    for i in patDF[xLabelCol]:
                        tmpdt = dt.datetime(2020, i, 1, 0, 0)
                        idxlab.append(tmpdt.strftime(xLabelFmt))
                        
                    
                    colr = self.Studies[scen].Color
                    lab = self.Studies[scen].ShortName
                    ls = self.Studies[scen].LineStyle
                    
                    if plotType=='avg':
                        ax.plot(idx, patDF.avg, lw=2.5, c=colr, label=lab, ls=ls)
                        ax.set_xticks(idx)
                        ax.set_xticklabels(idxlab)   
                        ax.yaxis.set_minor_locator(AutoMinorLocator())
                    
                    if plotType=='box':
                        dftmp['Month'] = [(m-9) if m>9 else m+3 for m in dftmp.index.month]
                        dftmp['Scen'] = [lab]*len(dftmp)
                        scen_df = pnd.concat([scen_df, dftmp])
                        pal[lab] = colr
                        labs.append(lab)
                        ax.yaxis.set_minor_locator(AutoMinorLocator())
                        
                    if plotType=='all':
                        data_cols = [y for y in patDF.columns if str(y).upper() not in ['CALYRMONTH','MONTH','YEAR','AVG','MAX','MIN']]
                        data_cols = [y for y in data_cols if 'q' not in str(y)]
                        quant_cols = [y for y in patDF.columns if 'q' in str(y)]
                        #s = set(quant_cols)
                        #data_cols2 = [x for x in data_cols if x not in s]
                        if 'lw' in kwargs:
                            lw = kwargs['lw']
                        else:
                            lw = 0.5
                        for c in data_cols:
                            if len(exclude_months)>0:
                                #tmpidx = pnd.Int64Index([i for i in idx if i not in  exclude_months])
                                patDFtmp = dcp(patDF[c])
                                patDFtmp[patDF.Month.isin(exclude_months)] = np.nan
                            else:
                                patDFtmp = dcp(patDF[c])
                                
                            if 'step' in kwargs:
                                if kwargs['step'] not in ['pre','post','mid']:
                                    swhere='pre'
                                else:
                                    swhere=kwargs['step']
                                ax.step(idx, patDFtmp, lw=lw, c=colr, alpha=0.2,
                                        where=swhere,label='tmp')   
                            else:
                                ax.plot(idx, patDFtmp, lw=lw, c=colr, alpha=0.2,
                                        label='tmp')

                            
                        hdl, labl = ax.get_legend_handles_labels()
                        print(hdl)
                        nhdls.append(hdl[0])

                        #nlabls.append('Position Analysis Realizations')
                        
                        if len(quant_cols)>0:
                            quantLW = 2.5
                            for c in quant_cols:
                                if len(exclude_months)>0:
                                    #tmpidx = pnd.Int64Index([i for i in idx if i not in  exclude_months])
                                    patDFtmp = dcp(patDF[c])
                                    patDFtmp[patDF.Month.isin(exclude_months)] = np.nan
                                else:
                                    patDFtmp = dcp(patDF[c])
                                
                                qlab = '%s pctile' %(float(c[1:])*10)
                                
                                if 'step' in kwargs:
                                    if kwargs['step'] not in ['pre','post','mid']:
                                        swhere='pre'
                                    else:
                                        swhere=kwargs['step']   
                                    ax.step(idx, patDFtmp, lw=quantLW, c=colr,
                                            alpha=1., label=qlab,
                                            where=swhere)
                                else:
                                    ax.plot(idx, patDFtmp, lw=quantLW, c=colr, alpha=1., label=qlab)
                                hdl, labl = ax.get_legend_handles_labels()
                                nhdls.append(hdl[-1])
                                nlabls.append(labl[-1])
                                quantLW = quantLW*1.8
                                
                        
                        ax.set_xticks(idx)
                        ax.set_xticklabels(idxlab)
                        
                        ax.yaxis.set_minor_locator(AutoMinorLocator())
                            
                    # the following code is for plotting all traces on the same
                    # monthly axis - not quite set up correctly for this
                    # TODO: fix all-trace plotting method
#                    for c in patDF.columns:
#                        if type(c)==str and c.upper() not in ['MONTH','YEAR', 'DAY']:
#                            ax.plot(idx, patDF.loc[:,c], color=colr,
#                                    label=lab, lw=1.2)
#                        elif type(c)!=str:
#                            ax.plot(idx, patDF.loc[:,c], color=colr,
#                                    label=lab, lw=1.2)
#                        else:
#                            pass
 
                    #return([fig,ax])
            if plotType=='box':
                g= sns.boxplot(x="Month", y=scen_df.columns[0],
                            hue="Scen", palette=pal,
                            data=scen_df,
                            saturation = 1,
                            fliersize=1,
                            linewidth=0.1,
                            medianprops={'color':'white','linewidth':1.0},
                            #linecolor='k',
                            #labels=idxlab,
                            #whiskerprops={'edgecolor':'k'},
                            ax=ax)
                ax.set_xticklabels(idxlab)
                plt.legend([], [], frameon=False) # gets rid of automatic seaborn legend
                
            legcols = nscens%2 + nscens//2
            
            ax.set_title(title, fontsize=16)
            ax.set_ylabel(ylabel + ' (%s)' %newUnits, fontsize=14)
            
            if nscens>5: #len(nhdls)>5:
                #plt.legend(nhdls, nlabls, bbox_to_anchor=(0.2,-0.25, 0.6, 0.15), ncol=legcols,frameon=False, fontsize=11)
                fig.legend(nhdls, nlabls, loc="outside lower center",
                           fontsize=11, ncol=legcols, frameon=False )
            else:
                #plt.legend(bbox_to_anchor=(0.2,-0.25, 0.6, 0.15), ncol=legcols,frameon=False, fontsize=11)
                fig.legend(loc="outside lower center",fontsize=13, frameon=False, ncol=nscens )
            
            plt.subplots_adjust(bottom=0.17,top=0.90) #,hspace=0.55)
            sns.despine()
            fig_dict[var] = [fig, ax]        
    
        return(fig_dict)
        
        
        
    
    def plot_multi_annualTS(self,title, annualize_on='A-Sep', how='auto-sum',
                            highlight_years=[],highlight_crit=True,**kwargs):
        '''
            make an annual time series plot with all scenarios in analysis for
            given variable var across multiple panels by 20-year periods
        
        '''
        fig_dict = {} # container to hold returned figure and associated data
        
        if 'exclude' in kwargs:
            exclude_scens = kwargs['exclude']
        else:
            exclude_scens = []
        
        hl_wy  = { '1920-1940': [],
                  '1940-1960': [],
                  '1960-1980': [],
                  '1980-2000': [],
                  '2000-2020': []}
        
        
        if len(highlight_years)>0:
            for yrs in  highlight_years:
                wy = yrs
                bdate = dt.date(wy-1,10,1)
                edate = dt.date(wy,10,1)
                if wy < 1940:
                    d = '1920-1940'
                elif (wy >= 1940) & (wy <1960):
                    d='1940-1960'
                elif (wy >=1960) & (wy < 1980):
                    d='1960-1980'
                elif (wy >=1980) & (wy < 2000):
                    d='1980-2000'
                else:
                    d='2000-2020'
                    
                hl_wy[d].append([wy, bdate, edate])
                
            
        elif highlight_crit:
            if self.Studies[0].WYT_FP=='DVfile':
                s = self.Studies[0]
                mask = s.WaterYearTypes==5
                tmp = s.WaterYearTypes.loc[mask.values, idx[:,'WYT_SAC_','WATERYEARTYPE']]
                tmp = tmp[tmp.index.month==9]
            else:
                if len(self.Studies)==1:  # if there's only one study, force it to be the baseline
                    s = self.Studies[0]
                    tmp = s.WaterYearTypes[s.WaterYearTypes.SACINDEX==5].WATERYEAR
                else:
                    for i,s in self.Studies.items():
                        if s.Baseline:
                            tmp = s.WaterYearTypes[s.WaterYearTypes.SACINDEX==5].WATERYEAR
                    
            if self.Studies[0].WYT_FP=='DVfile':     
                for i, val in tmp.iterrows():
                    wy = i.year
                    bdate = dt.date(wy-1,10,1)
                    edate = dt.date(wy,10,1)
                    if wy < 1940:
                        d = '1920-1940'
                    elif (wy >= 1940) & (wy <1960):
                        d='1940-1960'
                    elif (wy >=1960) & (wy < 1980):
                        d='1960-1980'
                    elif (wy >=1980) & (wy < 2000):
                        d='1980-2000'
                    else:
                        d='2000-2020'
                        
                    hl_wy[d].append([wy, bdate, edate])
            else:
                for crityrs in tmp:
                    wy = crityrs
                    bdate = dt.date(wy-1,10,1)
                    edate = dt.date(wy,10,1)
                    if wy < 1940:
                        d = '1920-1940'
                    elif (wy >= 1940) & (wy <1960):
                        d='1940-1960'
                    elif (wy >=1960) & (wy < 1980):
                        d='1960-1980'
                    elif (wy >=1980) & (wy < 2000):
                        d='1980-2000'
                    else:
                        d='2000-2020'
                        
                    hl_wy[d].append([wy, bdate, edate])
                
        else:
            hl_wy = hl_wy
            
        for var, dictDF in self.Staged.items():   # iterate through the staged variables, each one having a set of dataframes
                       
            with sns.plotting_context('paper', font_scale=1.3):
                fig, ax = plt.subplots(len(decades),1,sharex=False, sharey=True, figsize=(8,12))
        
                nscens = 0
                for scen, df in dictDF.items():   # iterate through the scenarios and dataframes
                    print(scen)
                    if scen in exclude_scens:
                        next
                    else:
                        nscens+=1
                        
                    origUnits = df.columns.get_level_values('Units')[0]
                    
                    data1 = df.copy()  
                     
                    scen_df = pnd.DataFrame()
    
                #for n, scen in self.Studies.items():
                    #data1 = self.Analytes[var][n] # data for analyte is indexed teh same as the studies
                    data2 = annualize(data1, on=annualize_on, how=how, colindex=[0])
                    
                    colr = self.Studies[scen].Color
                    lab = self.Studies[scen].ShortName
                    
                    if 'labels' in kwargs:
                        if kwargs['labels'].lower() in ['desc','description','descr']:
                            lab = self.Studies[scen].Description
                        else:
                            raise NotImplementedError('labels options limited to "description" for now')
                    
                    if self.Studies[scen].Baseline:
                        isBaseline=True
                        lw=2
                    else:
                        isBaseline=False
                        lw=1.2
                
                    for i,d in enumerate(decades):
                        data3 = data2.loc[decades[d][0]:decades[d][1]]
                        data = data3.shift(periods=-1, freq=annualize_on).append(data3.iloc[-1]) # shift the data so the step plot covers the water year
                         
                        # select water years for highlighting
                        hilite_yrs = hl_wy[d]
                        
                        for yr in hilite_yrs:
                            ax[i].axvspan(yr[1], yr[2], alpha=0.1, facecolor='0.8',edgecolor="None")
                        
                        ax[i].step(data.index, data.iloc[:,0], where='post', 
                          color=colr,
                          label=lab,lw=lw)
                        
                        ax[i].xaxis.set_minor_locator(mdates.YearLocator())
                        
                        if 'ytickfactors' in kwargs:
                            ax[i].yaxis.set_major_locator(MultipleLocator(kwargs['ytickfactors']))
                        ax[i].yaxis.set_minor_locator(AutoMinorLocator())
                        
                        ax[i].set_ylabel('Annual volume (TAF)', fontsize=10)
                        
                        # if d=='2000-2020':
                        #     ax[i].set_xlim((dt.date.toordinal(dt.date(1998,10,31)),dt.date.toordinal(dt.date(2020,9,30))))
                            
                        # if d=='1920-1940':
                        #     ax[i].set_xlim((dt.date.toordinal(dt.date(1918,10,31)),dt.date.toordinal(dt.date(1940,9,30))))
        
                        ax[i].set_xlabel('Year')
                        
                        if d=='2000-2020':
                         #   ax[i].set_xlim((dt.date.toordinal(dt.date(1998,10,31)),dt.date.toordinal(dt.date(2020,9,30))))
                             ax[i].set_xlim((dt.date(1998,10,31)),(dt.date(2020,9,30)))
                            
                        if d=='1920-1940':
                        #    ax[i].set_xlim((dt.date.toordinal(dt.date(1918,10,31)),dt.date.toordinal(dt.date(1940,9,30))))
                            ax[i].set_xlim((dt.date(1918,10,31)),(dt.date(1940,9,30)))
        
                    
            plt.suptitle(title, fontweight='bold',fontsize=20)
            
            #plt.legend(bbox_to_anchor=(0.2,-0.51, 0.6, 0.15), ncol=3,frameon=False, fontsize=11)
            axbox = ax[i].get_position()
            #+0.5*axbox.width
            # plt.legend(bbox_to_anchor=[axbox.x0-0.05*axbox.width, axbox.y0-0.06,0.6, 0.15], 
            #            bbox_transform=fig.transFigure,
            #            ncol=3,frameon=False, fontsize=10)
            legcols = nscens%2 + nscens//2
            h1, l1 = ax[0].get_legend_handles_labels()
            fig.legend(h1, l1, loc="outside lower center",fontsize=13, ncol=legcols )
            plt.subplots_adjust(bottom=0.15,top=0.95,hspace=0.55)
            sns.despine()
            
            fig_dict[var] = [fig, ax]
#            if SCRATCH:
#                plt.savefig(os.path.join(scratch_dir, 'climate_scenarios_ShastaInflowAnn_compare.png'), dpi=300)
        return(fig_dict)

    def plot_multi_monthlyTS(self,title,annualize_on='A-Sep',
                                highlight_years=[],highlight_crit=True,**kwargs):
            '''
                make an monthly time series plot with all scenarios in analysis for
                given variable var across multiple panels by 20-year periods
            
            '''
            if 'exclude' in kwargs:
                exclude_scens = kwargs['exclude']
            else:
                exclude_scens = []
            
            hl_wy  = { '1920-1940': [],
                      '1940-1960': [],
                      '1960-1980': [],
                      '1980-2000': [],
                      '2000-2020': []}
            
            # identify variable type for setting ylabel
            if 'vartype' in kwargs:
                vartype =kwargs['vartype']
            else:
                vartype = 'S' # assume storage if nothing else
            
            
            if len(highlight_years)>0:
                for yrs in  highlight_years:
                    wy = yrs
                    bdate = dt.date(wy-1,10,1)
                    edate = dt.date(wy,10,1)
                    if wy < 1940:
                        d = '1920-1940'
                    elif (wy >= 1940) & (wy <1960):
                        d='1940-1960'
                    elif (wy >=1960) & (wy < 1980):
                        d='1960-1980'
                    elif (wy >=1980) & (wy < 2000):
                        d='1980-2000'
                    else:
                        d='2000-2020'
                        
                    hl_wy[d].append([wy, bdate, edate])
                    
                
            elif highlight_crit:
                if self.Studies[0].WYT_FP=='DVfile':
                    s = self.Studies[0]
                    mask = s.WaterYearTypes==5
                    tmp = s.WaterYearTypes.loc[mask.values, idx[:,'WYT_SAC_','WATERYEARTYPE']]
                    tmp = tmp[tmp.index.month==9]
                else:
                    if len(self.Studies)==1:  # if there's only one study, force it to be the baseline
                        s = self.Studies[0]
                        tmp = s.WaterYearTypes[s.WaterYearTypes.SACINDEX==5].WATERYEAR
                    else:
                        for i,s in self.Studies.items():
                            if s.Baseline:
                                tmp = s.WaterYearTypes[s.WaterYearTypes.SACINDEX==5].WATERYEAR
                        
                if self.Studies[0].WYT_FP=='DVfile':     
                    for i, val in tmp.iterrows():
                        wy = i.year
                        bdate = dt.date(wy-1,10,1)
                        edate = dt.date(wy,10,1)
                        if wy < 1940:
                            d = '1920-1940'
                        elif (wy >= 1940) & (wy <1960):
                            d='1940-1960'
                        elif (wy >=1960) & (wy < 1980):
                            d='1960-1980'
                        elif (wy >=1980) & (wy < 2000):
                            d='1980-2000'
                        else:
                            d='2000-2020'
                            
                        hl_wy[d].append([wy, bdate, edate])
                else:
                    for crityrs in tmp:
                        wy = crityrs
                        bdate = dt.date(wy-1,10,1)
                        edate = dt.date(wy,10,1)
                        if wy < 1940:
                            d = '1920-1940'
                        elif (wy >= 1940) & (wy <1960):
                            d='1940-1960'
                        elif (wy >=1960) & (wy < 1980):
                            d='1960-1980'
                        elif (wy >=1980) & (wy < 2000):
                            d='1980-2000'
                        else:
                            d='2000-2020'
                            
                        hl_wy[d].append([wy, bdate, edate])
                    
            else:
                hl_wy = hl_wy
              
            fig_dict = {}
            for var, dictDF in self.Staged.items():   # iterate through the staged variables, each one having a set of dataframes
               
                 
                scen_df = pnd.DataFrame()
                pal = {}
                labs = []
                
                nlabls = []
                nhdls = []
            
                nscens = 0
                with sns.plotting_context('paper', font_scale=1.2):
                    fig, ax = plt.subplots(len(decades),1,sharex=False, sharey=True, figsize=(8,12))
            
                    for scen, df in dictDF.items():   # iterate through the scenarios and dataframes

                        if scen in exclude_scens:
                            continue
                        else:
                            nscens+=1

                        origUnits = df.columns.get_level_values('Units')[0]
                        
                        data1 = df.copy()  
                        
                        if vartype=='D':
                            if origUnits.upper() == 'CFS':
                                data1 = cfs_to_taf(data1)
    
                        if self.Studies[scen].Baseline:
                            isBaseline=True
                            lw=2
                        else:
                            isBaseline=False
                            if 'lw' in kwargs:
                                lw = kwargs['lw']
                            else:
                                lw=1.2
                    #for n, scen in self.Studies.items():
                        #data1 = self.Analytes[var][n] # data for analyte is indexed teh same as the studies
                        #data2 = annualize(data1, on=annualize_on, how=how, colindex=[0])
                        
                        #colr = scen.Color
                        #lab = scen.ShortName
 
                        if ('color' in kwargs) and (not isBaseline):
                            print("changing color")
                            colr = kwargs['color']
                        else:
                            colr = self.Studies[scen].Color
                            
                        lab = self.Studies[scen].ShortName
                        
                        if 'labels' in kwargs:
                            if kwargs['labels'].lower() in ['desc','description','descr']:
                                lab = self.Studies[scen].Description
                            else:
                                raise NotImplementedError('labels options limited to "description" for now')
                        
                        ls = self.Studies[scen].LineStyle
                        
                        maxv = np.nanmax(data1)
                        minv = np.nanmin(data1)
                        maxvlg = round(np.log10(maxv),0)
                        minvlg = np.log10(minv)
                        
                        for i,d in enumerate(decades):
                            data3 = data1.loc[decades[d][0]:decades[d][1]]
                            #data = data3.shift(periods=-1, freq=annualize_on).append(data3.iloc[-1]) # shift the data so the step plot covers the water year
                            data = data3 #.shift(periods=-1, freq=annualize_on).append(data3.iloc[-1]) # shift the data so the step plot covers the water year                            
                            
                            # select water years for highlighting
                            hilite_yrs = hl_wy[d]
                            
                            for yr in hilite_yrs:
                                ax[i].axvspan(yr[1], yr[2], alpha=0.1, 
                                              facecolor='0.8',edgecolor="None")
                            
                            if 'step' in kwargs:
                                ax[i].step(data.index, data.iloc[:,0], where='post', 
                                           color=colr,
                                           label=lab,lw=lw, ls=ls)
                            else:
                                ax[i].plot(data.index, data.iloc[:,0], 
                                           color=colr,
                                           label=lab,lw=lw, ls=ls)
                                
                            ax[i].xaxis.set_minor_locator(mdates.YearLocator())
                            
                            if 'ytickfactors' in kwargs:
                                ax[i].yaxis.set_major_locator(MultipleLocator(kwargs['ytickfactors']))
                            ax[i].yaxis.set_minor_locator(AutoMinorLocator())
                            
                            if vartype=='S':
                                ax[i].set_ylabel('Reservoir Storage (TAF)', fontsize=10)
                            elif vartype=='D':
                                ax[i].set_ylabel('Delivery (TAF)')
                            elif vartype=='CF':
                                ax[i].set_ylabel('Flow (CFS)')
                            elif vartype=='DI':
                                ax[i].set_ylabel('Diversion (CFS)')
                            else:
                                ax[i].set_ylabel('Volume')
                            
                            if d=='2000-2020':
                             #   ax[i].set_xlim((dt.date.toordinal(dt.date(1998,10,31)),dt.date.toordinal(dt.date(2020,9,30))))
                                 ax[i].set_xlim((dt.date(1998,10,31)),(dt.date(2020,9,30)))
                                
                            if d=='1920-1940':
                            #    ax[i].set_xlim((dt.date.toordinal(dt.date(1918,10,31)),dt.date.toordinal(dt.date(1940,9,30))))
                                ax[i].set_xlim((dt.date(1918,10,31)),(dt.date(1940,9,30)))
            
                            ax[i].set_xlabel('Year')
                            
                    plt.suptitle(title, fontweight='bold')
                    
                    #plt.legend(bbox_to_anchor=(0.2,-0.51, 0.6, 0.15), ncol=3,frameon=False, fontsize=11)
                    axbox = ax[i].get_position()
                    #+0.5*axbox.width
                    legcols = nscens%2 + nscens//2
                    h1, l1 = ax[0].get_legend_handles_labels()
                    # plt.legend(bbox_to_anchor=[axbox.x0-0.05*axbox.width, axbox.y0-0.06,0.6, 0.15], 
                    #            bbox_transform=fig.transFigure,
                    #            ncol=3,frameon=False, fontsize=10)
                    fig.legend(h1, l1, loc="outside lower center",fontsize=13,
                               ncol = legcols)
                    plt.subplots_adjust(bottom=0.15,top=0.95,hspace=0.55)
                    sns.despine()
                            
        #            if SCRATCH:
        #                plt.savefig(os.path.join(scratch_dir, 'climate_scenarios_ShastaInflowAnn_compare.png'), dpi=300)
                #return([fig,ax])        
                fig_dict[var] = [fig, ax]
            return(fig_dict)
    
    def plot_wyt_freq(self, wyt_var='SACINDEX',uses_cam_fcst=True):
        '''
            create grouped bar plot showing frequency of water year types
            
            options for wyt_var are:
                'SACindex', 'SJRindex', 'Shastaindex', 'AmerD893',
                'Featherindex', 'Trinityindex'
        '''
        

            
        title_dict = {'SACINDEX': 'Sacramento Water Year Types',
                      'SJRINDEX': 'San Joaquin Water Year Types',
                      'SHASTAINDEX': 'Shasta Water Year Index',
                      'AMERD893': 'American Water Year Index',
                      'FEATHERINDEX': 'Feather Water Year Index',
                      'TRINITYINDEX': 'Trinity Water Year Index'}
        
        fcst_vs_perffsight = {'SACINDEX': 'WYT_SAC_',
                              'SJRINDEX': 'WYT_SJR_',
                              'SHASTAINDEX':'WYT_SHASTA_CVP_',
                              'AMERD893': 'WYT_AMERD893_CVP_',
                              'FEATHERINDEX':'WYT_FEATHER_',
                              'TRINITYINDEX':'WYT_TRIN_'}
        title_dict_fcsts = {'WYT_SAC_':'Sacramento Water Year Types',
                            'WYT_SJR_':  'San Joaquin Water Year Types',
                            'WYT_SHASTA_CVP_': 'Shasta Water Year Index',
                            'WYT_AMERD893_CVP_': 'American Water Year Index',
                            'WYT_FEATHER_':'Feather Water Year Index',
                            'WYT_TRIN_': 'Trinity Water Year Index'}
       
        if uses_cam_fcst:
            wyt_var = fcst_vs_perffsight[wyt_var].upper()
            plt_title = title_dict_fcsts[wyt_var]
        else:
            wyt_var = wyt_var.upper()
            plt_title = title_dict[wyt_var]
        
        
        
        num_studies = len(self.Studies)
        tmp = {}
        for n, scen in self.Studies.items():
            if uses_cam_fcst:
                tmp[n]= [(len(scen.WaterYearTypes.loc[:,idx[:,wyt_var]].iloc[:,0].unique())),
                           sorted(scen.WaterYearTypes.loc[:,idx[:,wyt_var]].iloc[:,0].unique())]
            else:
                tmp[n]= [(len(scen.WaterYearTypes[wyt_var].unique())),
                       sorted(scen.WaterYearTypes[wyt_var].unique())]
        maxcat = 0
        maxcatID = 0
        for k in tmp:
            if tmp[k][0]>maxcat:
                maxcat=tmp[k][0]
                maxcatID = k
        num_cats = maxcat #max([i[0] for i in tmp.values()])
        cat_labels = tmp[maxcatID][1]
        print(num_cats)
        print(cat_labels)
        
        if wyt_var == 'SACINDEX' or wyt_var=='WYT_SAC_':
            cat_labels = ['Wet (1)','AN (2)','BN (3)','Dry (4)','Crit (5)']
        else:
            cat_labels = cat_labels # [scen.WaterYearTypes[wyt_var].unique()]
                
        
        pos= list(range(1,num_cats+1)) #[1, 2, 3, 4, 5]
        width= 1/(num_studies+1)
        with sns.plotting_context('paper', font_scale=1.2):
            fig, ax = plt.subplots(1,1, figsize=(10,5))
            
              #['Wet (1)','AN (2)','BN (3)','Dry (4)','Crit (5)']
            for n,i in enumerate(self.Studies):
                scen = self.Studies[n]
                

                print([p + width*n for p in pos])
                if uses_cam_fcst:
                    sel_idx = scen.WaterYearTypes.index.month==7
                    print(scen.WaterYearTypes.loc[sel_idx,idx[:,wyt_var]].iloc[:,0].value_counts(sort=False))
                    
                    ax.bar([p + width*n for p in pos],
                           scen.WaterYearTypes.loc[sel_idx,idx[:,wyt_var]].iloc[:,0].value_counts(sort=True),
                           width = width,
                           alpha=0.8,
                           color=scen.Color,
                           label=scen.ShortName)
                else:
                    print(scen.WaterYearTypes[wyt_var].value_counts(sort=False))
                    ax.bar([p + width*n for p in pos],
                           scen.WaterYearTypes[wyt_var].value_counts(sort=False),
                           width = width,
                           alpha=0.8,
                           color=scen.Color,
                           label=scen.ShortName)
                    
            ax.set_ylabel('Frequency of Year Type', fontsize=14)
            ax.set_title(plt_title, fontweight='bold', fontsize=20)
            
            ax.set_xticks([p +2.5*width for p in pos])
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            ax.set_xticklabels(cat_labels, fontsize=14)
            ax.yaxis.set_tick_params(labelsize=14)
            
            plt.legend(loc='upper right', ncol=2, frameon=False)
            sns.despine()
            
        return([fig,ax])
        
    # def plot_hist_(self, annualize_on=None):
    #     '''
    #         create histogram of data
            
    #         options for wyt_var are:
    #             'SACindex', 'SJRindex', 'Shastaindex', 'AmerD893',
    #             'Featherindex', 'Trinityindex'
    #     '''
        
        
class cs_study():

    def __init__(self, csObj,shortname, baseline=False,color='k', getWYT=True, **kwargs):
        self.CalSimObj = csObj
        self.Baseline = baseline
        self.Color= color
        self.ShortName = shortname
                
        if 'fullname' in kwargs:
            self.FullName = kwargs['fullname']
        else:
            self.FullName = ''
        
        if 'desc' in kwargs:
            self.Description = kwargs['desc']
        else:
            self.Description=''
        
        if 'linestyle' in kwargs:
            self.LineStyle = kwargs['linestyle']
        else:
            self.LineStyle = '-'
            
        if 'wytFP' in kwargs:
            if os.path.exists(kwargs['wytFP']):
                self.WYT_FP = kwargs['wytFP']
            else:
                print("specified water year type file doesn't exist at %s" %(kwargs['wytFP']))
                print("Trying the default wyt file path in the CalSim tables folder")
                
                
                wytFP2 = os.path.join(os.path.dirname(self.CalSimObj.LaunchFP), 'CONV','Run','Lookup','wytypes.table')
                if os.path.exists(wytFP2):
                    self.WYT_FP = wytFP2
                else:
                    print("Nope...that didnt' work either...going to have to try something different")
                    self.WYT_FP = None
        else:
            
            if csObj.Reorg:
                # check for CAM-calculated WYT in DV file
                if getWYT:
                    csObj.DVdata.getDVts(filter=['/WYT_SAC_/'])
                    self.WaterYearTypes = csObj.DVdata.DVtsDF
                    self.WYT_FP = 'DVfile' #os.path.join(os.path.dirname(self.CalSimObj.LaunchFP),'Run','Lookup','wytypes.table')
            else:
                wytFP2 = os.path.join(os.path.dirname(self.CalSimObj.LaunchFP), 'CONV','Run','Lookup','wytypes.table')
                if os.path.exists(wytFP2):
                    self.WYT_FP = wytFP2
                else:
                    print(f"-----Couldn't find WYT file at {wytFP2}\n\n")
                    self.WYT_FP = None
            
                if getWYT and self.WYT_FP != None:
                    tmpdf = cs3.read_CSlookup(self.WYT_FP)
                    colnames = [c.upper() for c in tmpdf.columns]
                    tmpdf.columns = colnames
                    self.WaterYearTypes = tmpdf #cs3.read_CSlookup(self.WYT_FP)
            
    def getSV(self, filter=''):
        if filter == '':
            self.CalSimObj.SVdata.getSVts()
        else:
            if type(filter)==list:
                self.CalSimObj.SVdata.getSVts(filter=filter)
            else:
                self.CalSimObj.SVdata.getSVts(filter=[filter])
        
    def getDV(self, filter=''):
        if filter =='':
            self.CalSimObj.DVdata.getDVts()
        else:
            if type(filter)==list:
                self.CalSimObj.DVdata.getDVts(filter=filter)
            else:
                self.CalSimObj.DVdata.getDVts(filter=[filter])
        
            
def cfs_to_taf(df):
    
    origcols = dict(zip(df.columns.names, df.columns[0]))
    cols = origcols.copy()
    cols['Units'] = 'TAF'
    coltup = tuple(cols.values())
    
    
    if type(df.index) in [pnd.core.indexes.datetimes.DatetimeIndex]:
        days_in_month = df.index.day
        fac = (86400./43560.)*days_in_month/1000.  
        
    df.loc[:,coltup] = df.loc[:,tuple(origcols.values())]*fac
    
    return(df.loc[:,[coltup]])

def taf_to_cfs(df):
    
    origcols = dict(zip(df.columns.names, df.columns[0]))
    cols = origcols.copy()
    cols['Units'] = 'CFS'
    coltup = tuple(cols.values())
    
    
    if type(df.index) in [pnd.core.indexes.datetimes.DatetimeIndex]:
        days_in_month = df.index.day
        fac = 1000/days_in_month * (43560.0/86400.0)  
        
    df.loc[:,coltup] = df.loc[:,tuple(origcols.values())]*fac
    
    return(df.loc[:,[coltup]])
    
def annualize(df, on='A-Sep', how='auto-sum', **kwargs):
    '''
        only works with a single column index/name specified for now
        TODO: deal with multiple columns to aggregate
    '''
    
    # if units info is provided in df, check this to inform 
    # the annual aggregation method
    units = df.columns.get_level_values('Units')[0]
    
    if units.upper() == 'CFS':
        if how.lower()=='sum':
            print("Warning! Units of data are rate (%s) but aggregation specified as sum" %units)
            print("Warning! Results will not be meaningful")
        if how.lower()[0:4]=='auto':
            df = cfs_to_taf(df)
            if len(how)>4:
                how_xtra = how.split('-')[1]
                how = how_xtra
            else:
                how='sum'
    else:
        if how.lower()[0:4]=='auto':
            if len(how)>4:
                how_xtra = how.split('-')[1]
                how = how_xtra
            else:
                how= 'sum'
            
    
    if 'colindex' in kwargs:
        howdict = {}
        cols = [df.columns[c] for c in kwargs['colindex']]
        for c in kwargs['colindex']:
            howdict[df.columns[c]] = how
            print(howdict)
            print(cols)
    elif 'colnames' in kwargs:
        howdict = {}
        cols = [c for c in kwargs['colnames']]
        for c in cols: #kwargs['colnames']:
            howdict[c] = how
    else:
        howdict={}
        cols = [df.columns[0]]
        if how !='sum':
            howdict[cols[0]] = how
        else:
            howdict[cols[0]] = 'sum'
    
    tmpdict = {}
    print(f'\nthis is the howdict: {howdict}\n')
    for c in cols:       
        print(f'Column: {c}')
        tmpdict[c] = df.loc[:,c].resample(on).apply(howdict[c]) #[c] #[c] #[c for c in cols]
    
    df_ann = pnd.concat(tmpdict, axis=1, names=('A','B','C','E','F', 'Type','Units')) #, names=[k for k in tmpdict.keys()])
    #df_ann = df[cols].resample(on).apply(howdict) #[c for c in cols]

    return(df_ann)

def single_exceed(df, col):
    series = df.loc[:, col]
    tmp = list(af.calc_exceed(series))
    probs = [i[0] for i in tmp]
    srtdvals = [i[1] for i in tmp]  
    name = series.name
    # initialize an empty dataframe based on time index and expected
    # DSS path name components
#    col_mindex = pnd.MultiIndex(levels=[[]]*7,
#                         codes=[[]]*7,
#                         names=[u'A', u'B', u'C',u'E',u'F',u'Type',u'Units'])
    newdf = pnd.DataFrame(index=probs, data=srtdvals, columns=[name])
    return(newdf)


def patternizer(df, freq='M', **kwargs):
    
    if type(df) ==pnd.Series:
        tt1 = pnd.DataFrame(data=df.copy(), columns=['value'])
    else:
        tt1 = df.copy()      

    
    
    if freq=='M':    
        tt1['CalYrMonth'] = tt1.index.month
        tt1['WYmonthorder'] = tt1.index.map(lambda x: af.wymo(x))
        tt1['Year'] = tt1.index.year
        tt1['WY'] = tt1.index.map(lambda x: af.addWY(x))
        
    if 'reference' in kwargs:
        if kwargs['reference']=='WY':    
            ttpv = pnd.pivot(tt1, index='WYmonthorder', columns='WY')[tt1.columns[0]] #['value']
            ttpv['Month'] = [10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        elif kwargs['reference']=='MarOct':  # for Sac temp mgmt season
            ttpv = pnd.pivot(tt1, index='CalYrMonth', columns='Year')[tt1.columns[0]]
            ttpv['Month'] = [1,2,3,4,5,6,7,8,9,10,11,12]
    else:
        # assume water year
        ttpv = pnd.pivot(tt1, index='WYmonthorder', columns='WY')[tt1.columns[0]] #['value']
        ttpv['Month'] = [10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        
    if 'summary' in kwargs:
        summary_meth = kwargs['summary']
        ret1 = ttpv.loc[:, ttpv.columns[0:-2]].agg(summary_meth, axis=1)
        ret2 = pnd.DataFrame(ret1)
        ret2['Month'] = ttpv['Month']
        return([ret2, ttpv])
    
    return(ttpv)        


def read_playbook(yamlFP):
    import yaml
    '''
        input: full filepath to valid yaml file specifying one or more sets of
               actions to take, with the adequate information to accompany eac
               action; 
               
               example: for copying records from one file to another, you'd include
                       the paths to the source and destination files along with the
                       list of record paths to copy
    '''
    plybk = []
    with open(yamlFP,'r') as of:
        tmp = yaml.safe_load_all(of)  #, Loader=yaml.SafeLoader)
        for t in tmp:
            plybk.append(t)
    return(plybk)

def run_playbook(plybk,base_dir='D:', **kwargs):
    
    # NOTE - if the paths in the playbook dictionary are relative, need
    # to provide the base directory that should preface all 
    # launch/config file paths
          
    if type(plybk)==dict:
        plybk = [plybk]
    else:
        plybk = plybk
       
    return_dict = Odict()
    print("made it here")
        
    if 'getWYT' in kwargs:
        getWYT = kwargs['getWYT']
    else:
        getWYT = True
        
    # get the setup and studies actions first
    ply_studies_keys = [i for i in range(len(plybk)) if plybk[i]['action_type'].lower() in ['setup', 'studies']]
    plt_plot_keys = [i for i in range(len(plybk)) if plybk[i]['action_type'].lower() in ['plot']]

#    for ply in plybk: # iterate through each of multiple actions if provided        
#        if ply['action_type'].lower() in ['setup', 'studies']:
    for psk in ply_studies_keys:

        #import cs3
        study_list = []
        ply = plybk[psk]
        for sty in ply['studies']:
            
            shtName = sty #list(sty.keys())[0]
            print("processing study ID: %s" %shtName)
            tmp = ply['studies'][sty] #sty[shtName]
            tkys = list(tmp.keys())
            
            if 'reorg' not in tmp.keys() or tmp['reorg']==True:
                reorg=True
            else:
                reorg=False
            
            if 'include' not in tmp.keys() or tmp['include']==True:

                if 'launchFP' in tmp:
                    launchFP = tmp['launchFP']
                    if not os.path.isabs(launchFP):
                        newlaunchFP = (Path(base_dir)/ launchFP).resolve() #os.path.join(base_dir, launchFP)).resolve()
                    else:
                        newlaunchFP = launchFP
                        
                    if not os.path.exists(newlaunchFP):
                        raise FileNotFoundError(f"Couldn't find launch file {newlaunchFP}")
                        return
                        
                    ldFLAG = False
                    if launchFP == None or launchFP=='':
                        ldFLAG = True
                        try:
                            launchDict = tmp['launchDict']
                            ldkys = launchDict.keys()
                            
                        except:
                            print("Insufficient launch information provided for study ID: %s" %shtName)
                            print("Need at least a launch file or a path to an SV file")
                            continue

                    cs2_cs3 = tmp['cs2cs3']
                    if ldFLAG:
                        tmpCalsim = cs3.calsim(launchDict=launchDict, 
                                              calsim_version=cs2_cs3, reorg=reorg)
                    else:
                        tmpCalsim = cs3.calsim(launchFP=newlaunchFP, reorg=reorg)
                elif 'configFP' in tmp:
                    cs2_cs3 = tmp['cs2cs3']
                    tmpCalsim = cs3.calsim(configFP = tmp['configFP'], reorg=reorg)
                    
                else:
                    print("No config or launch file provided - can't do anything")
                    return
                    
                baseline = tmp['baseline']
                color = tmp['color']
                
                desc = tmp['description'] 

                if 'ls' not in tmp:
                    linestyle = '-'
                else:
                    linestyle = tmp['ls']
                
                tmpStudy = cs_study(tmpCalsim, shtName, baseline=baseline, 
                                    color=color, desc = desc, linestyle=linestyle,
                                    getWYT= getWYT)
                study_list.append(tmpStudy)
                
            else:
                print("-->  skipping study ID: %s" %shtName)
            
        return_dict['analysis'] = cs_analysis(study_list)
            
    if len(plt_plot_keys) < 1:
        return(return_dict)
        
#    for plp in plt_plot_keys:
        
def plot_exceed(df, title, annualize_on='A-Sep', how='auto-sum',bl_id = "",
                           month_filter=[],annual_filter=[],
                           annual_filter_type="WY",
                           annual_filter_label="", 
                           return_data = False,
                           leg_fac1 = 0.2, # factor to center legend
                           **kwargs):
    '''
        month_filter: will plot exceedance based on just months meeting filter criteria (e.g. month_filter=[9] will include only September)
        annual filter: will plot exceedance based on just years meeting filter criteria
        
        
        kwargs:
            'reverse_x': reverses x-axis (exceedance probabilities) so that 100% is at the left
            'line_colors': dict containing scenario id (keys) and assigned color (values)
    '''
    fig_dict = {}
    
    if 'exclude' in kwargs:
        exclude_scen = kwargs['exclude']
    else:
        exclude_scen = []
        
    if type(title)!=list:
        title = [title]
        
    excd_data = {}
    var_cntr = 0
    # for var in df.columns:   # iterate through the scenarios in df columns
        
    var_data_ = {}
    var = np.unique([c[1] for c in df.columns])[0]
    with sns.plotting_context('notebook', font_scale=1.2):
        fig, ax = plt.subplots(1,1, figsize=(10,8))

        nscens = len(df.columns)
        colrs = sns.color_palette(palette='deep', n_colors=nscens, desat=0.9)
        for ix, scen in enumerate(df.columns): #, df in dictDF.items():   # iterate through the scenarios and dataframes
            scen_name = scen[0] #self.Studies[scen].ShortName
            
            if (scen in exclude_scen) or (scen_name in exclude_scen):
                continue
            
            if scen[0]==bl_id: #self.Studies[scen].Baseline:
                isBaseline=True
                lw=2
                thiszorder = 9999
                colr='k'
                if 'line_colors' in kwargs:
                    colr = kwargs['line_colors'][bl_id]
            else:
                isBaseline=False
                if 'lw' in kwargs:
                    lw = kwargs['lw']
                else:
                    lw=1.2
                thiszorder = 0
            
            ls = "-" #self.Studies[scen].LineStyle
            
            origUnits = df.columns.get_level_values('Units')[0]
            
            dftmp = df.loc[:, [scen]].copy()
            dftmp['WY'] = dftmp.index.map(af.addWY)
            
            #dftmp = dftmp.loc['1921-10-31':'2021-09-30',:]
            # FIRST: filter out jsut the years or months desired
            
            # check for annual filters - this filters out specific years
            if annual_filter == []:
                annfilt = dftmp
            else:
                if annual_filter_type.upper() in ["WY","WATERYEAR","WATER YEAR"]:
                    annfilt = dftmp[dftmp.WY.isin(annual_filter)]
                else:
                    annfilt = dftmp[dftmp.index.year.isin(annual_filter)]  #??? This only filters on calendar year - 
                                                                 # TODO: will need to add a way to do water year or delivery year
                                                                 # type filters
            
            # check for month filters - this selects for specific months
            if month_filter == []:
                monfilt = annfilt
            else:
                monfilt = annfilt[annfilt.index.month.isin(month_filter)]
                
            # drop the WY column to avoid the hassle of carrying that around
            #print(monfilt.head())
            filt = monfilt.drop(columns=['WY'])
            
            # SECOND: now do the aggregation/annualization
            if annualize_on != None:
                agg_df= annualize(filt, on=annualize_on, how=how, 
                                  colindex=[0], )
                finUnits = agg_df.columns.get_level_values('Units')[0]
            else:
                agg_df = filt
                finUnits = origUnits
#                        if annual_filter == []:
#                            ann_df= annualize(dftmp, on=annualize_on, how=how, colindex=[0])
#                        else:
#                            if annual_filter_type.upper() in ["WY","WATERYEAR","WATER YEAR"]:
#                                df2 = dftmp[dftmp.WY.isin(annual_filter)]
#                            else:
#                                df2 = dftmp[dftmp.index.year.isin(annual_filter)]  #??? This only filters on calendar year - 
#                                                                         # TODO: will need to add a way to do water year or delivery year
#                                                                         # type filters
                    #ann_df = annualize(df2, on=annualize_on, how=how, colindex=[0])
                                 
#                        if month_filter==[]:
#                            ann_df = annualize(df, on=annualize_on, how=how, colindex=[0])
#                        else:
#                            df2 = df[df.index.month.isin(month_filter)]
#                            ann_df = annualize(df2, on=annualize_on, how=how, colindex=[0])
                            
            excd_df = single_exceed(agg_df, agg_df.columns[0])
            
            var_data_[scen_name] = excd_df    
            
            if not isBaseline:
                if 'line_colors' in kwargs:
                    colr = kwargs['line_colors'][scen[0]]
                else:
                    colr = colrs[ix] #[self.Studies[scen].Color
                
            if 'labels' in kwargs:
                lab = kwargs['labels'][scen[0]]
            else:
                lab = scen[0] #self.Studies[scen].ShortName
            
            # if 'color' in kwargs:
            #     if type(kwargs['color'])==dict:
            #         if scen_name in kwargs['color']:  # by name
            #             colr = kwargs['color'][scen_name]
            #         elif scen in kwargs['color']: # by index
            #             colr = kwargs['color'][scen]
            #         if not isBaseline:
            #             lw = 1
            #     else:
            #         if not isBaseline:
            #             colr = kwargs['color']
                        
            if 'zorder' in kwargs:
                if type(kwargs['zorder'])==dict:
                    if scen_name in kwargs['zorder']:
                        thiszorder=kwargs['zorder'][scen_name]
                        
            ax.plot(excd_df.index, excd_df.iloc[:,0], color= colr, 
                    label=lab, lw=lw, zorder=thiszorder, ls=ls)
            
            ax.set_xlabel('Exceedance probability')
            
            if 'reverse_x' in kwargs:
                if kwargs['reverse_x']:
                    ax.invert_xaxis()
                    
            ax.xaxis.set_major_formatter(PercentFormatter(1.0))
            ax.xaxis.set_major_locator(MultipleLocator(0.1))
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            
            if finUnits.upper() in ['TAF','AF','ACRE-FEET']:
                ax.set_ylabel('Volume (%s)' %finUnits)
            else:
                ax.set_ylabel('Annual Average Flow (%s)' %finUnits)   #TODO: this is just a placeholder - need to add ability to deal with cfs, EC, other units and aggregations
            
            
            if finUnits.upper() =='TAF':
                ax2 = ax.secondary_yaxis( -0.15, functions=(taf_to_mcm, mcm_to_taf))
                ax2.yaxis.set_major_locator(AutoLocator()) #tkr.MultipleLocator(1))
                ax2.yaxis.set_minor_locator(AutoMinorLocator()) #tkr.MultipleLocator(0.2))
                ax2.set_ylabel('Million cubic meters',fontsize=14)
                plt.subplots_adjust(left=0.22)
            
            if month_filter != []:
                monthTitle = "months: "
                for m in month_filter:
                    if len(month_filter)==1:
                        monthTitle = f"{dt.date(2020,m,1):%B}"
                    else:
                        monthTitle += "%s, " %month_dict[m]
            else:
                monthTitle = ""
            
            if annual_filter != [] and annual_filter_label !="":
                annTitle = annual_filter_label
            elif annual_filter != []:
                annTitle = "Select years"
            else:
                annTitle = ""
                
            if annTitle=="":
                sep=''
            else:
                sep='-'
            if month_filter != []:
                full_title = title[var_cntr] + " Exceedance\n%s %s %s" %(annTitle,sep, monthTitle)
            else:
                full_title = title[var_cntr] + " Annual Exceedance\n%s %s %s" %(annTitle,sep, monthTitle)
            ax.set_title(full_title, fontsize=17, fontweight='bold')
    excd_data[var] = var_data_
    var_cntr+=1
        
    if nscens > 15:
        h1, l1 = ax.get_legend_handles_labels()
        #hcolors = list(set([hh.get_color() for hh in h1]))
        h2dict = {}
        for hh in h1:
            h2dict[hh.get_color()] = hh
        h2 = [v for k,v in h2dict.items()]
        if 'group_labels' not in kwargs:
            l2 = [l1[0]] + ['realizations']
        else:
            l2 = kwargs['group_labels']
        #l2 = []
        #h2 = [h1[0]]+[h1[1]] # get baseline + first realization handles
        #
        plt.legend(h2,l2, bbox_to_anchor=(0.2,-0.28, 0.6, 0.15), ncol=2,
                   frameon=False, fontsize=14)
    else:
        legcols = nscens%2 + nscens//2
        #plt.legend(bbox_to_anchor=(0.3,-0.25, 0.6, 0.15), ncol=legcols,frameon=False, fontsize=11)
        axbox = ax.get_position()
        #+0.5*axbox.width
        # plt.legend(bbox_to_anchor=[axbox.x0+leg_fac1*axbox.width, axbox.y0-0.14,0.6, 0.15], 
        #            bbox_transform=fig.transFigure,
        #            ncol=legcols,frameon=False, fontsize=13)
        fig.legend(loc='outside lower center', ncols=legcols, frameon=False, fontsize=13)
    plt.subplots_adjust(bottom=0.15)
    sns.despine()
    fig_dict[var] = [fig, ax]        
#            if SCRATCH:
#                plt.savefig(os.path.join(scratch_dir, 'climate_scenarios_ShastaInflowAnn_compare.png'), dpi=300)
    if return_data:
        return([fig_dict, excd_data])
    else:
        return(fig_dict)        
    
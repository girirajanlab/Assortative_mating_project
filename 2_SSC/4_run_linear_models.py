import pandas as pd
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

import scipy.stats as stats
import statsmodels.api as sm

matplotlib.rcParams['pdf.fonttype'] = 42

df=pd.read_csv('SSC_AM_family_table.csv')
df=df[df.linear_models=='X']

# Get input and response variables
covariates=['PC'+str(i) for i in range(1, 11)]+['Age']
inputs=['Sex', 'Autism PRS', 'Tier S SFARI SNV', 'SNV burden', 'Del (bp)', 'Dup (bp)', 'Biparental mean SRS', 'Biparental mean BAPQ']
phenos=['ADOS: social affect',
       'ADOS: restricted and repetitive behavior domain total score',
       'ADI: communication (verbal) domain score',
       'ADI: restricted and repetitive behavior domain total score',
       'ADI: social domain total score', 'RBS',
       'Parent-reported SRS: total raw scores', 'SCQ',
       'Insistence of sameness factor', 'Social interaction factor',
       'Sensory-motor behavior factor', 'Self-injurious behavior factor',
       'Idiosyncratic repetitive speech and behavior',
       'Communication skills factor',
       'Vineland Adaptive Behavior Scales: composite standard scores',
       'Full-scale IQ', 'Verbal IQ', 'Nonverbal IQ',
       'Developmental Coordination Disorders Questionnaire',
       'CBCL 2-5 internalizing', 'CBCL 2-5 externalizing',
       'CBCL 6-18 internalizing', 'CBCL 6-18 externalizing']

# Make sure all probands have input variables
df.dropna(subset=inputs+covariates, how='any', inplace=True, axis=0)

# Make sex and Tier S SFARI SNV binary
df['Sex']=df.Sex.map({'male':1, 'female':0})
df['Tier S SFARI SNV']=df['Tier S SFARI SNV'].map({'LOF':1, '0.0':0, 'Missense':1, 'Splice':1})

def plot_scatter(subdf, xlab, ylab, ax, norm=False):
    x=subdf[xlab].to_list()
    y=subdf[ylab].to_list()
    if xlab in ['Sex', 'Tier S SFARI mutation']:
        sns.violinplot(data=subdf, x=v, y=ylab, ax=ax)
        r, p = stats.ttest_ind(x, y)
        test='T'
        y=0.1
    else:
        sns.scatterplot(data=subdf, x=xlab, y=ylab, ax=ax, alpha=0.2)
        r, p = stats.pearsonr(x, y)
        test='R'
        y=0.5
    if p < 0.01:
        ax.annotate("%s: %.2f\np: %.3E" % (test, r, p), xy = (0.1, y), xycoords = ax.transAxes, bbox=dict(boxstyle="round", fc="0.8"))
    elif p > 0.05:
        ax.annotate("%s: %.2f\np: %.2f" % (test, r, p), xy = (0.1, y), xycoords = ax.transAxes, bbox=dict(boxstyle="round", fc="0.8"))
    else:
        ax.annotate("%s: %.2f\np: %.4f" % (test, r, p), xy = (0.1, y), xycoords = ax.transAxes, bbox=dict(boxstyle="round", fc="0.8"))
    ax.set_ylabel(ylab)
    ax.set_xlabel(xlab)
    if norm:
        ax.set_ylabel(ylab+'\n(normalized)')

def plot_res(res, names, pdf, title, covariates=None):
    fig, ax=plt.subplots(figsize=(3, 4))
    
    est=list(res.params)[1:]
    err=list(res.bse)[1:]
    p=list(res.pvalues)[1:]
    ci=res.conf_int(alpha=0.05)
    ci_lo=list(ci[0])[1:]
    ci_up=list(ci[1])[1:]
    r2=res.rsquared_adj
    
    names2=names.copy()
    est2=est.copy()
    err2=err.copy()
    ci_lo2=ci_lo.copy()
    ci_up2=ci_up.copy()
    p2=p.copy()
    
    if covariates!=None:
        for cov in covariates:
            idx=names2.index(cov)
            
            del names2[idx]
            del est2[idx]
            del err2[idx]
            del ci_lo2[idx]
            del ci_up2[idx]
            del p2[idx]
            
    ys=[i for i in range(len(names2))]
    ys.reverse()
    
    # Estimates and error bars
    # Mark nominally significant results with blue color
    plt.scatter([est2[i] for i in range(len(p2)) if p2[i]>0.05], [ys[i] for i in range(len(p2)) if p2[i]>0.05], color='k')
    plt.scatter([est2[i] for i in range(len(p2)) if p2[i]<=0.05], [ys[i] for i in range(len(p2)) if p2[i]<=0.05], color='#00AEEF')
    for i in range(len(ys)):
        color='k'
        if p2[i]<=0.05:
            color='#00AEEF'
        plt.plot([ci_lo2[i], ci_up2[i]], [ys[i], ys[i]], color=color)
    
    # Axis labels
    plt.yticks(ys, names2)
    
    # Line at 0
    plt.plot([0, 0], [-0.5, max(ys)+0.5], color='k', ls=':', zorder=0)
    plt.ylim(-0.25, max(ys)+0.25)
    
    # Add title
    plt.title(title)
    
    plt.tight_layout()
    pdf.savefig()
    plt.close()
    
    # Also parse stats and return
    out_stat=[]
    for i in range(len(ys)+len(covariates)):
        out_stat.append([names[i], est[i], err[i], ci_lo[i], ci_up[i], p[i], r2])
    
    return out_stat

def make_model(y, x):
    X=sm.add_constant(x)
    mod=sm.OLS(y, X)
    return mod.fit()

# Create a model for each phenotype
pdf=PdfPages('Figures/4_quant_models.pdf')
pdf2=PdfPages('Figures/4_quant_corr.pdf')
stat_lst=[]
for p2 in phenos:
    print(p2)
    subdf=df[(~df[p2].isnull())][[p2]+covariates+inputs].copy()
    subdf.dropna(how='any', inplace=True, axis=0)
    # Scale all inputs
    # Plot correlations before and after normalization
    fig, axs=plt.subplots(nrows=2, ncols=len(covariates+inputs), figsize=(22, 5))
    for i, v in enumerate(covariates+inputs):
        plot_scatter(subdf, v, p2, axs[0, i])
    subdf[p2]=(subdf[p2]-subdf[p2].mean())/subdf[p2].std()
    for i, v in enumerate(covariates+inputs):
        if v not in ['Sex', 'Tier S SFARI SNV']:
            if '(bp)' in v:
                subdf[v]=np.log10(df[v]+1)
            subdf[v]=(subdf[v]-subdf[v].mean())/subdf[v].std()
        plot_scatter(subdf, v, p2, axs[1, i])
    pdf2.savefig()
    plt.close()
    
    # Run linear model
    y=subdf[p2].to_numpy()
    x=subdf[covariates+inputs]
    mod=make_model(y, x)
    stat=plot_res(mod, covariates+inputs, pdf, p2, covariates=covariates)
    for s in stat:
        stat_lst.append([p2]+s+[subdf.shape[0]])

pdf.close()
pdf2.close()

stat_df=pd.DataFrame(stat_lst, columns=['proband_phenotype', 'variable', 'estimate', 'error', 'ci_lower', 'ci_upper', 'p', 'adj_r2', 'n'])
stat_df['FDR']=stat_df.p*stat_df.shape[0]
stat_df.loc[stat_df.FDR>1, 'FDR']=1

# Save to file
stat_df.to_csv('Result_tables/4_model_output.csv', index=False)
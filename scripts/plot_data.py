"""Render the retained quantitative figure families from exact CSV data."""
from fractions import Fraction as F
from pathlib import Path
import csv

def plot(data,figures):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import PercentFormatter
    def rows(name):
        with (data/(name+'.csv')).open() as f:return list(csv.DictReader(f))
    plt.rcParams.update({'font.size':10,'axes.labelsize':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'figure.facecolor':'white'})
    fig,ax=plt.subplots(figsize=(6.3,3.25),layout='constrained');rr=rows('cycle_noise');c=9/16
    for L,col,ls in zip((3,5,9,21),('#14688f','#008979','#ac711b','#945573'),('-','--','-.',':')):
        a=[r for r in rr if int(r['L'])==L];u=[float(F(r['u_exact'])) for r in a];y=[float(F(r['F_exact'])) for r in a];cut=(1-(2*c**(1/3)-1)**(1/L))/2
        ax.plot(u,y,color=col,ls=ls,lw=2,label=f'$L={L}$');ax.plot([cut],[c],'o',ms=4,color=col)
    ax.axhline(c,color='#50565c',lw=1,ls=(0,(4,4)),label='Free overlap $9/16$');ax.set(xlim=(0,.09),ylim=(.1,1.03),xlabel='Independent physical phase-flip probability',ylabel='Squared target fidelity');ax.xaxis.set_major_formatter(PercentFormatter(1,decimals=0));ax.grid(axis='y',alpha=.2);ax.legend(frameon=False,ncol=2,loc='upper right',fontsize=9);fig.savefig(figures/'cycle_noise.pdf');plt.close(fig)
    fig,ax=plt.subplots(figsize=(6.3,3.25),layout='constrained');rr=rows('batch_yield')
    for e,sty in zip((0,1,3),('o-','s--','^-.')):
        a=[r for r in rr if int(r['ell_representative'])==e];ax.plot([int(r['N']) for r in a],[int(r['yield']) for r in a],sty,label='1 or 2' if e==1 else str(e))
    ax.set(xlabel='Aligned source copies N',ylabel='Exact neutral CCZ yield',xticks=list(range(1,13)));ax.legend(title='Reference dimension',frameon=False);fig.savefig(figures/'batch_yield.pdf');plt.close(fig)

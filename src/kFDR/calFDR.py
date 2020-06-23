from statsmodels.stats import multitest
import pandas as pd

#---------------------------------------------------------------------------#
# https://www.statsmodels.org/dev/_modules/statsmodels/stats/multitest.html#multipletests
# https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html
#---------------------------------------------------------------------------#

def fdr_values(p_value):
    _b=multitest.multipletests(p_value, alpha=0.05, method='b', is_sorted=False, returnsorted=False)
    _b_v=pd.DataFrame(_b[1])
    _s=multitest.multipletests(p_value, alpha=0.05, method='s', is_sorted=False, returnsorted=False)
    _s_v=pd.DataFrame(_s[1])
    _h=multitest.multipletests(p_value, alpha=0.05, method='h', is_sorted=False, returnsorted=False)
    _h_v=pd.DataFrame(_h[1])
    _hs=multitest.multipletests(p_value, alpha=0.05, method='hs', is_sorted=False, returnsorted=False)
    _hs_v=pd.DataFrame(_hs[1])
    _sh=multitest.multipletests(p_value, alpha=0.05, method='sh', is_sorted=False, returnsorted=False)
    _sh_v=pd.DataFrame(_sh[1])
    _ho=multitest.multipletests(p_value, alpha=0.05, method='ho', is_sorted=False, returnsorted=False)
    _ho_v=pd.DataFrame(_ho[1])
    _fdr_bh=multitest.multipletests(p_value, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    _fdr_bh_v=pd.DataFrame(_fdr_bh[1])
    _fdr_by=multitest.multipletests(p_value, alpha=0.05, method='fdr_by', is_sorted=False, returnsorted=False)
    _fdr_by_v=pd.DataFrame(_fdr_by[1])
    _fdr_tsbh=multitest.multipletests(p_value, alpha=0.05, method='fdr_tsbh', is_sorted=False, returnsorted=False)
    _fdr_tsbh_v=pd.DataFrame(_fdr_tsbh[1])
    _fdr_tsbky=multitest.multipletests(p_value, alpha=0.05, method='fdr_tsbky', is_sorted=False, returnsorted=False)
    _fdr_tsbky_v=pd.DataFrame(_fdr_tsbky[1])
    _fdr_gbs=multitest.multipletests(p_value, alpha=0.05, method='fdr_gbs', is_sorted=False, returnsorted=False)
    _fdr_gbs_v=pd.DataFrame(_fdr_gbs[1])
    fdr_all= pd.concat([_b_v, _s_v, _h_v, _hs_v, _sh_v, _ho_v, _fdr_bh_v, _fdr_by_v, _fdr_tsbh_v, _fdr_tsbky_v, _fdr_gbs_v], axis=1, sort=False)
    fdr_all.columns = ['Bonferroni', 'Sidak', 'Holm', 'Holm-Sidak', 'Simes-Hochberg', 'Hommel', 'FDR Benjamini-Hochberg', 'FDR Benjamini-Yekutieli', 'FDR 2-stage Benjamini-Hochberg', 'FDR 2-stage Benjamini-Krieger-Yekutieli', 'FDR adaptive Gavrilov-Benjamini-Sarkar']
    return fdr_all



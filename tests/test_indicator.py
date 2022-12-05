from endorse.common import File, workdir,EndorseCache
from endorse.indicator import indicators
from endorse import plots


def test_indicator():
    EndorseCache.instance().expire_all()
    #pvd_file = File("test_data/trans_m_01/solute_fields.pvd")
    case = "trans_m_00"
    pvd_file = File(f"test_data/{case}/solute_fields.pvd")
    with workdir('sandbox'):
        inds = indicators(pvd_file, 'U235_conc', (-10, 10))
        plots.plot_indicators(inds, file=case)
        ind_time_max = [ind.time_max()[1] for ind in inds]
        print(ind_time_max)
import time
import numpy as np
#from mlmc.mc_level import Level


class MLMC:
    """
    Multilevel Monte Carlo method
    """

    def __init__(self, level_data, level_times):
        """
        :param levels: List[np.array], shape (n_level_samples, 2, n_variables)
        """
        self.level_data = level_data
        self.level_times = level_times

    @property
    def n_levels(self):
        """
        Number of levels
        """
        return len(self.levels)

    @property
    def n_samples(self):
        """
        Level samples
        """
        return np.array([len(l) for l in self.level_data])

    def level_variances(self):
        pass



    # @property
    # def n_nan_samples(self):
    #     """
    #     Level nan samples
    #     """
    #     return np.array([len(l.nan_samples) for l in self.levels])

    # @property
    # def sim_steps(self):
    #     return np.array([Simulation.log_interpolation(self.step_range, lvl.step) for lvl in self.levels])
    #
    # def sample_range(self, n0, nL):
    #     """
    #     Geometric sequence of L elements decreasing from n0 to nL.
    #     Useful to set number of samples explicitly.
    #     :param n0: int
    #     :param nL: int
    #     :return: np.array of length L = n_levels.
    #     """
    #     return np.round(np.exp2(np.linspace(np.log2(n0), np.log2(nL), self.n_levels))).astype(int)

    #
    # def subsample(self, sub_samples=None):
    #     """
    #     :param sub_samples: None - use all generated samples
    #                 array of ints, shape = n_levels; choose given number of sub samples from computed samples
    #     :return: None
    #     """
    #     if sub_samples is None:
    #         sub_samples = [None] * self.n_levels
    #     assert len(sub_samples) == self.n_levels, "{} != {}".format(len(sub_samples), self.n_levels)
    #     for ns, level in zip(sub_samples, self.levels):
    #         level.subsample(ns)
    #
    # def update_moments(self, moments_fn):
    #     for level in self.levels:
    #         level.evaluate_moments(moments_fn, force=True)
    #
    # def clean_levels(self):
    #     """
    #     Reset all levels
    #     :return: None
    #     """
    #     for level in self.levels:
    #         level.reset()
    #
    # def clean_subsamples(self):
    #     """
    #     Clean level subsamples
    #     :return: None
    #     """
    #     for level in self.levels:
    #         level.subsample(None)
    #
    # def get_sample_times(self):
    #     """
    #     The total average duration of one sample per level (fine + coarse together)
    #     :return: list
    #     """
    #     return [level.sample_time() for level in self.levels]

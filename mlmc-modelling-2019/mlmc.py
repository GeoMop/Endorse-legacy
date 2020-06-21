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
        self.level_data = [np.array(l) for l in level_data]
        self.level_times = np.array(level_times)

        # checks
        assert self.n_levels > 0
        assert self.level_times.shape == (self.n_levels, )
        for il, l in enumerate(self.level_data):
            assert l.shape[1] == 2
            if il > 0:
                assert l.shape[-1] == self.level_data[0].shape[-1]
    @property
    def n_variables(self):
        return self.level_data[0].shape[-1]

    @property
    def n_levels(self):
        """
        Number of levels
        """
        return len(self.level_data)

    @property
    def n_samples(self):
        """
        Level samples
        """
        return np.array([len(l) for l in self.level_data])

    def level_diff_variances(self):
        vars = []
        for l in self.level_data:
            # l shape (N,2,M)
            var_diff = np.var(l[:, 0, :] - l[:, 1, :], axis=0, ddof=1)
            vars.append(var_diff)
        return np.array(vars)    # shape (L, M)

    def level_variances(self):
        vars = []
        for l in self.level_data:
            # l shape (N,2,M)
            vars.append(np.var(l, axis=0, ddof=1))
        return np.array(vars)    # shape (L, 2, M)

    def level_means(self):
        """
        :return: (L, 2, M)
        """
        means = []
        for l in self.level_data:
            # l shape (N,2,M)
            means.append(np.mean(l, axis=0))
        return np.array(means)

    def variable_means(self):
        sums = []
        for l in self.level_data:
            # l shape (N,2,M)
            variable_sums = np.sum(l[:,0,:], axis=0)
            sums.append(variable_sums)
        n_fine_samples = np.sum(self.n_samples)
        return np.sum(sums, axis=0) / n_fine_samples


    def estimate_n_samples_for_target_variance(self, target_variance, prescribe_vars=None):
        """
        Estimate optimal number of samples for individual levels that should provide a target variance of
        resulting moment estimate. Number of samples are directly set to levels.
        This also set given moment functions to be used for further estimates if not specified otherwise.
        TODO: separate target_variance per moment
        :param target_variance: Constrain to achieve this variance.
        :param prescribe_vars: vars[ L, M] for all levels L, and M variables can be weighted from `level_variances` using `variable_means`.
        :return: np.array with number of optimal samples for individual levels, shape (L,)
        """
        vars = prescribe_vars
        assert vars.shape == (self.n_levels, self.n_variables)
        time_est = np.array(self.level_times)
        sqrt_var_n = np.sqrt(vars * time_est[:, None])  # moments in rows, levels in cols
        total = np.sum(sqrt_var_n, axis=0)  # sum over levels
        n_samples_estimate = np.round( sqrt_var_n / time_est[:, None] * total[None, :] / target_variance).astype(int)  # moments in cols
        n_samples_estimate_safe = np.minimum(n_samples_estimate, vars * self.n_levels / target_variance)
        n_samples_estimate_safe = np.maximum(n_samples_estimate_safe, 2)    # at least 2 samples for every level
        n_samples = np.max(n_samples_estimate_safe, axis=1).astype(int)     # maximum over variables
        return n_samples

    def estimate(self, n_samples, seed=0):
        """
        Subsample collected samples, using given `n_samples` vector.
        Perform MLMC estimate of variables using the subsample.
        :return: estimate_of_var_means, estimate_of_variance_of_estimate ; shapes (M,), (M,)
        """
        np.random.seed(seed)
        means = []
        vars = []
        for nl, l in zip(n_samples, self.level_data):
            if seed is None:
                subset = np.arange(nl)
            else:
                subset = np.random.choice(len(l), nl)
            l_sub = l[subset, :, :]
            l_diff = l_sub[:, 0, :] - l_sub[:, 1, :]
            l_mean = np.mean(l_diff, axis=0)
            l_var = np.var(l_diff, axis=0, ddof=1)
            means.append(l_mean)
            vars.append(l_var)
        means = np.sum(np.array(means), axis=0)
        vars = np.sum(np.array(vars) / n_samples[:, None], axis=0)
        return np.array(means), np.array(vars)



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

import openmc

class LeaprRun:
    """
    Converts vibrational densities of state for materials
    into S(alpha, beta) tables without the pain of using NJOY.

    You may include secondary scatterers, but only when using the free gas
    or diffusion approximations for their density of state. The short collision
    time treatment for secondary scatterers is currently not implemented, but
    could be with a little additional work.

    For an example use of this class, see the included Jupyter notebook
    in this repo.
    """

    _SECONDARY_SCATTERER_TREATMENTS = {'short-collision-time': 0,
                                       'free': 1,
                                       'diffusion': 2}

    _COHERENT_ELASTIC_CRYSTALS = {'graphite': 1,
                                  'beryllium': 2,
                                  'beryllium oxide': 3,
                                  'aluminum': 4,
                                  'lead': 5,
                                  'iron': 6}

    _SAB_TYPES = {'symmetric': 0,
                  'asymmetric': 1}

    # Currently unused, untested. Not hard to fill this out later.
    # _COLD_HYDROGEN_TYPES = {'ortho hydrogen':1,
    #                         'para hydrogen': 2,
    #                         'ortho deuterium': 3,
    #                         'para deuterium': 4}
    # _COLD_HYDROGEN_TREATMENT = {'vineyard': 1,
    #                             'skold': 2}

    def __init__(self):
        self.title = 'OpenMC LEAPR run'

        # must be added through add_temperature
        self._temperatures = []
        self._rho = []
        self._tbetas = []  # integral of DOS. computed automatically.
        self._delta_rho = None
        self._n_rho = None

        self._material_number = None
        self._zaid = None
        self._awr = None
        self._free_atom_xs = None
        self._alphas = None
        self._betas = None
        self._n_scattering_atoms = None
        self._coherent_elastic_crystal = None

        # whether a secondary scatterer is present, e.g. O in H2O.
        self._secondary_scatter_A = None
        self._secondary_scatter_free_xs = None
        self._num_secondary_scatterer_atoms = None
        self._secondary_scatter_treatment = None

        # Settings for Egelstaff-Schofield parameters used by liquids.
        # These are different at each temperature.
        self._twt = []
        self._egelstaff_schofield_c = []

        self._sab_type = 'symmetric'
        self._n_phonon = 200
        self._sout = 0

        # todo: add support for liquid hydrogen
        # self._cold_hydrogen_type = None
        # self._cold_hydrogen_treatment = None

    @property
    def material_number(self):
        return self._material_number

    @material_number.setter
    def material_number(self, matno):
        cv.check_type('material number', matno, int)
        self._material_number = matno

    @property
    def zaid(self):
        return self._zaid

    @zaid.setter
    def zaid(self, x):
        assert x > 0
        self._zaid = x

    @property
    def awr(self):
        return self._awr

    @awr.setter
    def awr(self, x):
        cv.check_type('atomic weight ratio', x, float)
        assert x > 0
        self._awr = x

    @property
    def free_atom_xs(self):
        return self._free_atom_xs

    @free_atom_xs.setter
    def free_atom_xs(self, xs):
        cv.check_type('free atom cross section', xs, float)
        assert xs > 0
        self._free_atom_xs = xs

    @property
    def alphas(self):
        return self._alphas

    @alphas.setter
    def alphas(self, x):
        try:
            self._alphas = np.array(x, dtype=float)
        except ValueError:
            raise ValueError("Cannot convert alpha grid to float.")

    @property
    def betas(self):
        return self._betas

    @betas.setter
    def betas(self, x):
        try:
            self._betas = np.array(x, dtype=float)
        except ValueError:
            raise ValueError("Cannot convert betas grid to float.")

    @property
    def n_scattering_atoms(self):
        return self._n_scattering_atoms

    @n_scattering_atoms.setter
    def n_scattering_atoms(self, n):
        cv.check_type('number of scattering atoms', n, int)
        assert n > 0
        self._n_scattering_atoms = n

    @property
    def temperatures(self) -> list[float]:
        return self._temperatures

    @temperatures.setter
    def temperatures(self, Tvals):
        raise Exception("Set temperatures through the add_temperature method.")

    @property
    def secondary_scatterer_A(self):
        return self._secondary_scatter_A

    @secondary_scatterer_A.setter
    def secondary_scatterer_A(self, A):
        cv.check_type('secondary scatterer AWR', A, float)
        assert A > 0.0
        self._secondary_scatter_A = A

    @property
    def secondary_scatterer_free_xs(self):
        return self._secondary_scatter_free_xs

    @secondary_scatterer_free_xs.setter
    def secondary_scatterer_free_xs(self, xs):
        cv.check_type('secondary scatterer free xs', xs, float)
        assert xs > 0.0
        self._secondary_scatter_free_xs = xs

    @property
    def num_secondary_scatterer_atoms(self):
        return self._num_secondary_scatterer_atoms

    @num_secondary_scatterer_atoms.setter
    def num_secondary_scatterer_atoms(self, n):
        cv.check_type('number secondary scatterer atoms', n, int)

        # keep as none if zero passed in
        if n == 0:
            return
        assert n > 0
        self._num_secondary_scatterer_atoms = n

    @property
    def twt(self):
        return self._twt

    @twt.setter
    def twt(self, w):
        raise Exception(
            "must set translational weight through add_temperature")

    @property
    def egelstaff_schofield_c(self):
        return self._egelstaff_schofield_c

    @egelstaff_schofield_c.setter
    def egelstaff_schofield_c(self, c):
        raise Exception(
            "must set egelstaff_schofield_c through add_temperature")

    @property
    def sab_type(self):
        return self._sab_type

    @sab_type.setter
    def sab_type(self, new_type):
        cv.checkvalue('S(a, b) type', new_type,
                      self.__class__._SAB_TYPES.keys())
        self._sab_type = new_type

    @property
    def n_phonon(self):
        return self._n_phonon

    @n_phonon.setter
    def n_phonon(self, n):
        cv.check_type('number of phonons', n, int)
        assert n > 0
        self._n_phonon = n

    @property
    def sout(self):
        return self._n_phonon

    @n_phonon.setter
    def sout(self, n):
        cv.check_type('mix choice of sab', n, int)
        assert n == 0 or n==1 or n==2
        self._sout = n

    @property
    def secondary_scatterer_treatment(self) -> str:
        return self._secondary_scatter_treatment

    @secondary_scatterer_treatment.setter
    def secondary_scatterer_treatment(self, treatment_type):
        if treatment_type == 'short-collision-time':
            raise NotImplemented("Have not implemented capability to do short collision "
                                 "time treatment of secondary scatterers. The LeaprRun class"
                                 "needs to be modified to write densities of state for the secondary "
                                 "scatterer at each temperature point.")
        cv.checkvalue('secondary scatterer treatment', treatment_type,
                      self.__class__._SECONDARY_SCATTERER_TREATMENTS)
        self._secondary_scatter_treatment = treatment_type

    def set_rho_grid(self, rho_grid) -> None:
        rho_diff = np.diff(rho_grid)
        if not np.all(rho_diff == rho_diff[0]):
            raise Exception("Density of states grid must be uniform.")
        else:
            self._delta_rho = rho_diff[0]
            self._n_rho = len(rho_grid)

    @property
    def n_rho(self):
        return self._n_rho

    @n_rho.setter
    def n_rho(self, n):
        cv.check_type('n_rho', n, int)
        assert n > 0
        self._n_rho = n

    @property
    def delta_rho(self):
        return self._delta_rho

    @delta_rho.setter
    def delta_rho(self, d):
        cv.check_type('delta_rho', d, float)
        assert d > 0
        self._delta_rho = d

    def add_temperature_point(self, T, rho_vals, twt=0.0, c=0.0) -> None:
        """
        Adds a density of states for temperature T.
        rho_vals are on the rho energy grid, in units of eV.
        """
        if self._delta_rho is None or self._n_rho is None:
            raise Exception("Must set rho grid before adding temperature points."
                            "Either set both n_rho and delta_rho or pass a rho grid"
                            " to set_rho_grid.")

        self._temperatures.append(T)

        # convert to numpy array if necessary
        if not isinstance(rho_vals, np.ndarray):
            rho_vals = np.array(rho_vals)

        self._rho.append(rho_vals)
        self._tbetas.append(np.trapz(rho_vals, dx=self._delta_rho))

        assert twt >= 0.0
        assert c >= 0.0
        self._egelstaff_schofield_c.append(c)
        self._twt.append(twt)

    # returns the number corresponding to the coherent elastic crystal type
    def _elastic_crystal_code(self) -> int:
        if self._coherent_elastic_crystal is None:
            return 0
        else:
            return self.__class__._COHERENT_ELASTIC_CRYSTALS[self._coherent_elastic_crystal]

    def _has_all_secondary_scatterer_options(self) -> bool:
        secondary_scatter_options = (self._secondary_scatter_A,
                                     self.secondary_scatterer_free_xs,
                                     self._num_secondary_scatterer_atoms,
                                     self._secondary_scatter_treatment)

        # if any secondary scatterer option is present, all must be present.
        if any([v is not None for v in secondary_scatter_options]):
            if None in secondary_scatter_options:
                raise Exception("Not all secondary scatterer options were set."
                                "Must set all if one has been set.")
            return True
        else:
            return False

    def _addcard(self, *values):
        str_values = [str(v) for v in values]
        self._commands += ' '.join(str_values)
        self._commands += '/\n'

    def _make_commands(self) -> str:
        """
        Writes the NJOY input file to a string, which is returned
        from this function. Other stuff handles writing it to a file
        and running NJOY on it.
        """

        # internally for storing input deck
        self._commands = "leapr\n"

        self._addcard(24)  # use tape 24 (arbitrary)

        # card 2
        self._addcard("'", self.title, "'")

        # card 3
        if not self._temperatures:
            raise ValueError("No temperature points were added to the run.")
        self._addcard(len(self._temperatures), 1, self._n_phonon)

        # card 4
        if self._material_number is None:
            raise Exception("must set material number")
        if self._zaid is None:
            raise Exception("must set primary scatterer ZAID")
        self._addcard(self._material_number, self._zaid,
                      self.__class__._SAB_TYPES[self._sab_type],
                      0, 1e-75, self._sout)

        # card 5 (note: liquid hydrogen options to be added here)
        if self._awr is None:
            raise Exception("Must set primary scatterer ZAID")
        if self._free_atom_xs is None:
            raise Exception(
                "Must set primary scatterer free atom cross section")
        if self._n_scattering_atoms is None:
            raise Exception("Must set number of primary scatterer atoms")
        self._addcard(self._awr, self._free_atom_xs, self._n_scattering_atoms,
                      self._elastic_crystal_code())

        # card 6
        if self._has_all_secondary_scatterer_options():
            self._addcard(1, self.__class__._SECONDARY_SCATTERER_TREATMENTS[_self.secondary_scatterer_treatment],
                          self._secondary_scatter_A, self._secondary_scatter_free_xs,
                          second._num_secondary_scatterer_atoms)
        else:
            self._addcard(0)

        # card 7, 8, 9
        if self._alphas is None or self._betas is None:
            raise Exception("must set beta and alpha grid")
        self._addcard(len(self._alphas), len(self._betas))
        mlw = 70  # max line width
        self._addcard(np.array2string(self._alphas, max_line_width=mlw)[1:-1])
        self._addcard(np.array2string(self._betas, max_line_width=mlw)[1:-1])

        # The remaining cards are repeated for every temperature value
        for i_T, T in enumerate(self._temperatures):
            self._addcard(T)  # card 10
            self._addcard(self._delta_rho, self._n_rho)  # card 11
            self._addcard(np.array2string(
                self._rho[i_T], max_line_width=mlw)[1:-1])  # card 12
            self._addcard(
                self._twt[i_T], self._egelstaff_schofield_c[i_T], self._tbetas[i_T])  # card 13
            self._addcard(0)  # card 14 (no discrete oscillators)
        self._addcard("'made in openmc'")
        self._addcard("")  # blank
        self._commands += " stop\n"
        return self._commands

    def run(self, stdout=False, njoy_exec='njoy', output_filename='leapr_out'):
        """
        Runs LEAPR and outputs to an ENDF S(a, b) file called leapr_out.
        """
        commands = self._make_commands()
        openmc.data.njoy.run(commands, {}, {24: output_filename},
            stdout=stdout, njoy_exec=njoy_exec)

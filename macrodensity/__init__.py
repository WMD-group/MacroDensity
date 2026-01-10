"""
MacroDensity is a package to read, process and plot electrostatic potential and electron density
files from electronic structure calculations.
"""

from macrodensity.averages import (
    macroscopic_average,
    planar_average,
    spherical_average,
    travelling_volume_average,
    volume_average,
)
from macrodensity.density import element_vol, gradient_magnitude
from macrodensity.io import (
    get_band_extrema,
    read_gulp_potential,
    read_vasp_density,
    read_vasp_density_classic,
    read_vasp_parchg,
)
from macrodensity.plotting import (
    energy_band_alignment_diagram,
    plot_active_plane,
    plot_active_space,
    plot_field_at_point,
    plot_on_site_potential,
    plot_planar_average,
    plot_plane_field,
    plot_variation_along_vector,
)
from macrodensity.tools import (
    bulk_interstitial_alignment,
    bulk_vac,
    create_plotting_mesh,
    diff_potentials,
    dipole_correction,
    extend_potential,
    match_resolution,
    matched_spline_generate,
    scissors_shift,
    sort_potential,
    spline_generate,
    subs_potentials,
    translate_grid,
)
from macrodensity.utils import (
    GCD,
    GCD_List,
    density_2_grid,
    get_third_coordinate,
    get_volume,
    inverse_participation_ratio,
    matrix_2_abc,
    number_in_field,
    numbers_2_grid,
    one_2_2d,
    points_2_plane,
    vector_2_abscissa,
)

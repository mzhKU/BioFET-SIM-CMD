from bio_lib import *

# ************************************************************************
# CONDUCTIVITY
# The BioFET-SIM Compute Unit.
# ........................................................................  
def compute(rho, nw_len, nw_rad,
                 lay_ox, L_d, L_tf, lay_bf,
                 eps_1, eps_2, eps_3, n_0,
                 nw_type, num_prot):

    # Computing sensitivity of wire.
    gamma = Gamma(nw_rad, lay_ox, L_d, L_tf, eps_1, eps_2, eps_3) 
    
    # Multiple Charge Model critical feature.
    GammaSigma = 0.0
    # Here we run over the charge distribution.
    # mod.rho[i][2]: z-coordinate of q_i ([Ang]).
    # mod.rho[i][3]: q_i ([Elementary charge]).
    # lay_ox:        Oxide layer thickness ([nm]).
    # lay_bf:        Biofunctionalization layer thickness ([nm]).
    # li_tot:        Distance of charge site in protein to nano wire ([nm]).
    # offset:        Coordinate of lowest charge*-1.0 ([nm]).

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # <<UNIT>> 28.12.2011: Switching to [nm].
    #offset = get_z_offset(mod.rho)
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    offset = get_z_offset(rho)*0.1
    for i in range(len(rho)): 

        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # <<UNIT>> 28.12.2011: Switching to [nm].
        #li_tot      = lay_bf + lay_ox + float(mod.rho[i][2]) + offset# + pos
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        li_tot      = lay_bf + lay_ox + float(rho[i][2])*0.1 + offset# + pos
        gamma_li    = Gamma_li(nw_rad, li_tot, L_d) 
        sigma_bi    = Sigma_bi(nw_rad, li_tot, nw_len, num_prot, rho[i][3]) 
        GammaSigma += gamma_li*sigma_bi
    
    # Converts to meter because [q_elem] is Coulomb and [n_0] is in SI units.
    nw_rad = nw_rad*1E-9 

    # Change in relative conductivity.
    dG_G0 = 2.0/(nw_rad*q_elem*n_0)*gamma*GammaSigma
    
    if nw_type == 'P':
        dG_G0 *= -1 
    
    return dG_G0
# ------------------------------------------------------------------------ 

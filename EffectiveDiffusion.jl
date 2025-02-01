using LinearAlgebra

function General_Deff(YY, Temp, radius, poro, P)
    #=
    YY=[1 2 3 4 5 6 7 8 9 10];
    radius=5;
    poro=0.5;
    P=13;
    Temp=22;
    =#
   # Constants
    H2_index = 1;
    CO_index = 2;
    CH4_index = 3;
    CO2_index = 4;
    PH3_index = 5;
    H2O_index = 6;
    N2_index = 7;
    O2_index = 8;
    Ar_index = 9;
    He_index = 10;
    Mw = [2.0, 28.0, 16.0, 44.0, 34.0, 18.0, 28.0, 32.0, 40.0, 2.0]  # Molecular weights (kg/kmol)
    n_species = length(Mw)
    small = 1.0e-30
    great = 1.0e+30
    pi = π
    R = 8.314  # Universal gas constant (J/mol-K)
    convf = 1.0e3  # Conversion factor

    # Shape factor for tortuosity
    alpha = 0.75
    tao = (1 - alpha * (1 - poro))^(-1)
    sphericity = 1.0
    k = sphericity * (2 * radius)^2 / 150.0 * poro^3 / (1 - poro)^2
    sigma = zeros(n_species)
    epsik = zeros(n_species)
    # Lennard-Jones parameters (sigma: Ångstroms, epsik: K)
    #sigma = [2.827, 3.690, 3.758, 3.941, 2.641, 3.798, 3.981, 3.467, 3.35, 2.57]
    #epsik = [59.7, 91.7, 148.6, 195.2, 809.1, 71.4, 251.5, 106.7, 143.2, 10.8]
	sigma[H2_index] = 2.827;
	sigma[CO_index] = 3.690;
	sigma[CH4_index] = 3.758; 
	sigma[CO2_index] = 3.941;
	sigma[H2O_index] = 2.641; 
	sigma[N2_index] = 3.798;
	sigma[PH3_index] = 3.981;
    sigma[O2_index] = 3.467;
    sigma[Ar_index] = 3.35;
    sigma[He_index] = 2.57;

    epsik[H2_index] = 59.7;
	epsik[CO_index] = 91.7;
	epsik[CH4_index] = 148.6;
	epsik[CO2_index] = 195.2; 
	epsik[H2O_index] = 809.1;
	epsik[N2_index]= 71.4; 
	epsik[PH3_index] = 251.5;
    epsik[O2_index]= 106.7;
    epsik[Ar_index] = 143.2;
    epsik[He_index] = 10.8;
    # Chapman-Enskog parameters
    A = 1.06036; B = 0.15610; C = 0.19300; D = 0.47635
    E = 1.03587; F = 1.52996; G = 1.76474; H = 3.89411

    # Initialize outputs
    Dk = zeros(n_species)        # Knudsen diffusivity
    D_bin = zeros(n_species, n_species)  # Binary diffusivity
    Dm = zeros(n_species)        # Molecular diffusivity
    Deff = zeros(n_species)      # Effective diffusivity
    alpha_factor = zeros(n_species)

    # Calculate Knudsen diffusivity
    for i in 1:n_species
        Dk[i] = (2 / 3) * sqrt(8 * R * Temp * convf / (pi * Mw[i])) * radius
    end
    sigma_bin=zeros(n_species, n_species);
    epsik_bin=zeros(n_species, n_species);
    Tn=zeros(n_species, n_species);
    Omega_D=zeros(n_species, n_species);
    # Calculate binary diffusivity
    for i in 1:n_species
        for j in 1:n_species
            sigma_bin[i, j] = (sigma[i] + sigma[j]) / 2
            epsik_bin[i, j] = sqrt(epsik[i] * epsik[j])
            Tn[i, j] = Temp / (epsik_bin[i,j] + small)
            Omega_D[i, j] = (A / (Tn[i,j] + small)^B) + (C / exp(D * Tn[i,j] + small)) + (E / exp(F * Tn[i,j] + small)) + (G / exp(H * Tn[i,j] + small))
            D_bin[i, j] = (0.001858 / 10000.0) * sqrt((Temp^3 * (Mw[i] + Mw[j])) / (Mw[i] * Mw[j] + small)) /
                          (P * sigma_bin[i,j]^2 * Omega_D[i,j] + small)
        end
    end
    return sigma_bin,epsik_bin,Tn,D_bin;
    # Calculate molecular diffusion coefficient
    for i in 1:n_species
        sum_term = 0.0
        for j in 1:n_species
            if i != j
                sum_term += YY[j] / (D_bin[i, j] + small)
            end
        end
        Dm[i] = (1.0 - YY[i]) / (sum_term + small)
    end 
    Mm=0.0;
    for i in 1:n_species
        Mm=Mm + YY[i]*Mw[i];
    end
    # Calculate effective diffusivity

    for i in 1:n_species
        alpha_factor[i] = 1.0 - sqrt(Mw[i] / (Mm + small))
        Deff[i] = (poro / tao) * 1.0 / ((1.0 - alpha_factor[i] * YY[i]) / (Dm[i] + small) +1.0 / (Dk[i] + small) + small)
    end
    return Deff, alpha_factor
    return D_bin, Dm, Dk, tao, k
    
end

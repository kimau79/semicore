########################################### Information ############################################
# Tested on Julia Version 1.4.2

# Packages reqiured:
using SpecialFunctions 		#error function required for star profiles
using DelimitedFiles		#Optional, but more convenient to read and write arrays

# date last updated: 21 May 2021
# Command to run the code from Unix:
# julia main.jl (value_of_vk) (value_of_tau)
# example: "julia main.jl 20 3"

####################################################################################################

include("myFunctions_hb.jl")

############################################ Parameters ############################################
################### Constants ####################
const Mo = 1.98847e30  # Solar mass [kg]
const kpc = 3.08567758128e19  # Kiloparsec [m]
const h = 0.6727  # "little h" [dimensionless]. Defined by H = 100 * h [km s-1 Mpc-1]
const H = 0.1  # Hubble constant at present [h km s-1 kpc-1]
const G = 43007.1  # Gravitational constant [1e-10 kpc Mo-1 (km s-1)-2]
const rho_c = 3 * H ^ 2 / (8 * pi * G)  # Critical density [1e10 Mo kpc-3 h2]

############## NFW halo parameters ###############
const M_vir = 5.17  # [1e10 Mo h-1]
const c = 25.5 / 1.17  # Concentration parameter of the NFW halo
const rho_avg = 200 * rho_c  # Average density within virial radius [1e10 Mo kpc-3 h2]. rho_avg = delta_vir * rho_c

################# DDM parameters #################
const v_k = parse(Float64,ARGS[1])  # Recoil velocity of newborn daughter particles [km s-1]
const tau = parse(Float64,ARGS[2])  # Half-life of mother particles [Gyr]

############## Time step parameters ##############
const t_end = 14  # Ending time [Gyr]
const numOfSteps = 14  # Total number of time evolution intervals
const Additional_steps = 2  #If the final time step is too large, can further break into small steps

################ Shell parameters ################
const firstShellThickness = 1e-6  # Thickness of first shell [kpc h-1]. If interested range of radius starts from 1e-n, use 1e-(n-2) for good accuracy
const shellThicknessFactor = 1.0032  # Thickness of each consecutive shell grows exponentially by this rate factor. Originally 1.0032
const extend_factor = 4  # Maximum halo radius at initialization is set as R_vir * extend_factor. 4 is recommended

############### The Ratio of Stars in NFW Profile ########
Sprofile_type=2 		#0 = off; 1 = gaussian; 2 = singular isothermal sphere
star_to_DM_ratio = 1/100	#Usuing DM's virial mass (not total mass) as base, this ratio gives the mass of stars with the "virial" radius specified below
star_vir_radius = 0.1 		#Together with the parameter above, it specifies the mass of star desired within a radius. (in kpc)
radius_to_sd = 2		#(Gaussian) This number sets the ratio between the radius and SD of the gaussian distribution. It means that the SD = star_vir_radius / radius_to_sd
cutoff_SingIso = 0.5		#(Isothermal) since singular isothermal profile diverges, cutoff req'd (hard cutoff)


################# Miscellaneous ##################
const tol_ellipseGuess = 0.00001 / 100  # Tolerance for bisection method in ellipseRadii() [in the unit of shellThickness]. Smaller tolerance gives more accurate r_max and r_min. 0.00001 / 100 is recommended
const orderOfpolynomial = 14  # Order of polynomial for fitting res(x). 14 is recommended
####################################################################################################

######################################## Calculations ##############################################
# NFW profile parameters
rho_0 = rho_avg * c ^ 3 / 3 / (log(1 + c) - c / (1 + c))
R_vir = (3 * M_vir / (4 * pi * rho_avg)) ^ (1 / 3)
R_s = R_vir / c
#rho_0 = 1.6022248e-3 / h / h
#R_vir = 42.63597 * h
#R_s = 3 * h
#const c = R_vir / R_s
#const M_vir = 4 * pi * rho_avg / 3 * R_vir ^ 3
const NFW_params = [rho_0, R_s, c]
star_params=[star_to_DM_ratio, star_vir_radius, radius_to_sd]

# For time step array
delta_m = (1 - exp(-log(2) * t_end / tau)) / numOfSteps  # For having equal amounts of mass decayed every time step
t = zeros(numOfSteps + 1)  # Array containing time steps [Gyr]
t[1] = 0
for i in 2:numOfSteps + 1
    t[i] = - log(1 - (i - 1) * delta_m) * tau / log(2)
end

if Additional_steps!=0
t_temp = zeros(numOfSteps + Additional_steps + 1)
for i = 1:numOfSteps
t_temp[i]=t[i]
end
filled = numOfSteps
global numOfSteps+=Additional_steps
for j=1:Additional_steps+1
t_temp[filled + j]=t[filled]+(t[filled+1]-t[filled])/(Additional_steps+1)*j
end
global t = t_temp
end

# Total number of shells
numOfShells = floor(Int, log(1 - NFW_params[2] * NFW_params[3] * extend_factor  / firstShellThickness * (1 - shellThicknessFactor)) / log(shellThicknessFactor)) + 1
####################################################################################################

########################################## Algorithm ###############################################

    # For calculating runtime
    functionStart = time_ns()
    stepStart = time_ns()

    # Master folder
    folderName = "dmOnly"
    if !isdir(folderName)
        mkdir(folderName)
    end

    # Parameter-specific subfolder
    folderName = folderName * "/" * string(M_vir) * "_" * string(v_k) * "_" * string(tau) * "_" * string(numOfShells) * "_" * string(numOfSteps) * "_" * string(star_params)
run_id = -1
for i=1:100
if !isdir(folderName*"run"*string(i))
	global run_id = i
	break
end
end
folderName *="run"*string(run_id)
        mkdir(folderName)

    # Results subsubfolder
    folderName_results = folderName * "/" * "results"
    if !isdir(folderName_results)
        mkdir(folderName_results)
    end


    # Print ALL parameters to a file
    paramsFileName = folderName * "/params.txt"
    f = open(paramsFileName, "w")
    println(f, "Mo=", Mo)
    println(f, "kpc=", kpc)
    println(f, "h=", h)
    println(f, "H=", H)
    println(f, "G=", G)
    println(f, "rho_c=", rho_c)
    println(f, "M_vir=", M_vir)
    println(f, "c=", c)
    println(f, "rho_avg=", rho_avg)
    println(f, "v_k=", v_k)
    println(f, "tau=", tau)
    println(f, "t_end=", t_end)
    println(f, "numOfSteps=", numOfSteps)
    println(f, "firstShellThickness=", firstShellThickness)
    println(f, "shellThicknessFactor=", shellThicknessFactor)
    println(f, "extend_factor=", extend_factor)
    println(f, "tol_ellipseGuess=", tol_ellipseGuess)
    println(f, "orderOfpolynomial=", orderOfpolynomial)
    println(f, "rho_0=", rho_0)
    println(f, "R_s=", R_s)
    println(f, "R_vir=", R_vir)
    println(f, "NFW_params=", NFW_params)
    println(f, "delta_m=", delta_m)
    println(f, "t=", t)
    println(f, "numOfShells=", numOfShells)

println(t)


function dmOnly()
    ############# Initializing NFW halo ##############
    t_0 = t[1]
    println("Initializing at t=$t_0 Gyr...")

    # Display the main parameters
    println("tau=", tau,", v_k=", v_k,", c=", c)
    println("R_vir=", R_vir, ", R_s=", R_s, ", rho_0=", rho_0)
    println("extend_factor=", extend_factor, ", orderOfpolynomial=", orderOfpolynomial)

    # Initialize NFW total shells
    Tshells_radii, Tshells_mass = NFW_shells(NFW_params, numOfShells, shellThicknessFactor, extend_factor)
    # Initialize daughter shells (empty)
    Dshells_mass = zeros(size(Tshells_radii, 1))  # No daughters at initialization

    Sshells_mass = zeros(size(Tshells_radii, 1))

    # Initialize Star shells
	if Sprofile_type==1 #Gaussian
		star_cum_mass = M_vir * star_to_DM_ratio
		#star_cum_mass = 4*pi*rho_0*R_s^3*(log(1+star_vir_radius/R_s)+R_s/(R_s+star_vir_radius)-1)*star_to_DM_ratio
		k_star = 2*(star_vir_radius/radius_to_sd)^2
		A_star = star_cum_mass/2/sqrt(pi*k_star)/(sqrt(pi*k_star)/2*erf(star_vir_radius/sqrt(k_star))-star_vir_radius*exp(-star_vir_radius^2/k_star))
		Cum_mass_star1(r) = 2*A_star*sqrt(pi*k_star)*(sqrt(pi*k_star)/2*erf(r/sqrt(k_star))-r*exp(-r^2/k_star))
		Sshells_mass[1]=Cum_mass_star1(Tshells_radii[1,2])
		Cum_mass_temp = Sshells_mass[1]
		for i=2:size(Tshells_radii, 1)
			Sshells_mass[i]=Cum_mass_star1(Tshells_radii[i,2])-Cum_mass_temp
			Cum_mass_temp += Sshells_mass[i]
		end
	elseif Sprofile_type==2 #Isothermal
		star_cum_mass = M_vir * star_to_DM_ratio
		sigma_v_sq = star_cum_mass * G / star_vir_radius	#rms velocity of stars
		Cum_mass_star2(r) = sigma_v_sq / G * r
		Sshells_mass[1] = Cum_mass_star2(Tshells_radii[1,2])
		Cum_mass_temp = Sshells_mass[1]
		for i=2:size(Tshells_radii,1)
			if Tshells_radii[i,2]<=cutoff_SingIso
				Sshells_mass[i] = Cum_mass_star2(Tshells_radii[i,2])-Cum_mass_temp
				Cum_mass_temp += Sshells_mass[i]
			elseif Tshells_radii[i,1]<=cutoff_SingIso
				Sshells_mass[i] = Cum_mass_star2(cutoff_SingIso)
			else
				Sshells_mass[i]=0
			end
		end
	end



    # Initialize Mother shells
    Mshells_mass = Tshells_mass # At initialization all DM particles
    Tshells_mass += Sshells_mass

    Tshells_enclosedMass = enclosedMass(Tshells_radii, Tshells_mass)  # Enclosed mass of all particles
    Tshells_GPE = NFW_GPE(Tshells_radii, NFW_params, G)  # NFW GPE for this initial distribution

    # Intermediate results at each time step are outputted
    MfileName = folderName * "/M_t=$t_0.txt"
    printToFile(Tshells_radii, Mshells_mass, MfileName, G)
    DfileName = folderName * "/D_t=$t_0.txt"
    printToFile(Tshells_radii, Dshells_mass, DfileName, G)
    SfileName = folderName * "/S_t=$t_0.txt"
    printToFile(Tshells_radii, Sshells_mass, SfileName, G)
    TfileName = folderName * "/T_t=$t_0.txt"
    printToFile(Tshells_radii, Tshells_mass, TfileName, G)
    GPEfileName = folderName * "/GPE_t=$t_0"*".txt"
    printToFile_GPE(Tshells_radii, Tshells_GPE, GPEfileName)
    DMTfileName = folderName * "/DMT_t=$t_0.txt"
    printToFile(Tshells_radii, Mshells_mass+Dshells_mass, DMTfileName, G)

    stepResultsFileName = folderName * "/stepResults.txt"
    g = open(stepResultsFileName, "w")

    timeTaken = (time_ns() - stepStart) / 1e9
    println("Time taken for this step: ", timeTaken, "s\n")
    println(g, t_0, "\t", timeTaken)

    ############# Time evolution starts ##############
    for t_i in 2:numOfSteps + 1
        stepStart = time_ns()

        t_rounded = round(t[t_i], digits = 2)  # Rounding off t for readability in terminal
        println("Working on ", t_i - 1, " / $numOfSteps step (", t_rounded, " Gyr)...")  # Tracking progress

        # Proportion of undecayed mother particles
        p_undecayed = exp(log(1 / 2) * t[t_i] / tau) / exp(log(1 / 2) * t[t_i - 1] / tau)

        # Calculate (per mass) L and total E of mothers from the total mass distribution
        Mshells_L = L(Tshells_radii, Tshells_enclosedMass, G)
        Mshells_totalE_afterDecay = totalE_afterDecay(Tshells_radii, Tshells_GPE, Mshells_L, v_k)

        # Solve equation for r_max, r_min
        Mshells_ellipseRadii = ellipseRadii(Mshells_L, Mshells_totalE_afterDecay, Tshells_radii, Tshells_GPE, tol_ellipseGuess)

        # Decay the mothers in the shells and distribute the new daughters
        Mshells_mass, Dshells_decayedMass = updateShellsMass(Tshells_radii, Mshells_ellipseRadii, Mshells_mass, p_undecayed, Mshells_L, Mshells_totalE_afterDecay, Tshells_GPE, Tshells_enclosedMass, t_i, orderOfpolynomial, G, NFW_params)

        # Now we have:
        # 1. Mshells (remaining mothers)
        # 2. Dshells (daughters accumulated from before)
        # 3. Dshells_decayed (new daughters born in this time step)
	# 4. Sshells (not affected by decays)

        # Prepare total enclosed mass values for adiabatic expansion
        Dshells_mass = Dshells_mass + Dshells_decayedMass  # Combine newborn daughters with old daughters
        Tshells_mass_updated = Mshells_mass + Dshells_mass + Sshells_mass  # Total distribution
        Tshells_enclosedMass_updated = enclosedMass(Tshells_radii, Tshells_mass_updated)  # Total enclosed mass
        Tshells_GPE_updated = GPE(Tshells_radii, Tshells_mass_updated, Tshells_enclosedMass_updated, G)  # GPE for this total distribution (not useful)

        # Results before adiabatic expansion at each time step are outputted
        MfileName = folderName * "/M_beforeAdia_t=$t_rounded.txt"
        printToFile(Tshells_radii, Mshells_mass, MfileName, G)
        DfileName = folderName * "/D_beforeAdia_t=$t_rounded.txt"
        printToFile(Tshells_radii, Dshells_mass, DfileName, G)
        DdefileName = folderName * "/Dde_beforeAdia_t=$t_rounded.txt"
        printToFile(Tshells_radii, Dshells_decayedMass, DdefileName, G)
        SfileName = folderName * "/S_beforeAdia_t=$t_rounded.txt"
        printToFile(Tshells_radii, Sshells_mass, SfileName, G)
        TfileName = folderName * "/T_beforeAdia_t=$t_rounded.txt"
        printToFile(Tshells_radii, Tshells_mass_updated, TfileName, G)
        GPEfileName = folderName * "/GPE_beforeAdia_t=$t_rounded.txt"
        printToFile_GPE(Tshells_radii, Tshells_GPE_updated, GPEfileName)

        # Adiabatic expansions (applied to both mothers and daughters)
        Mshells_mass = adiabaticExpansion(Tshells_radii, Mshells_mass, Tshells_enclosedMass, Tshells_enclosedMass_updated)
        Dshells_mass = adiabaticExpansion(Tshells_radii, Dshells_mass, Tshells_enclosedMass, Tshells_enclosedMass_updated)
        Sshells_mass = adiabaticExpansion(Tshells_radii, Sshells_mass, Tshells_enclosedMass, Tshells_enclosedMass_updated)

        # Update total mass shells
        Tshells_mass = Mshells_mass + Dshells_mass + Sshells_mass
        Tshells_enclosedMass = enclosedMass(Tshells_radii, Tshells_mass)
        Tshells_GPE = GPE(Tshells_radii, Tshells_mass, Tshells_enclosedMass, G)

        # Results at the end of each time step are outputted
        MfileName = folderName * "/M_t=$t_rounded.txt"
        printToFile(Tshells_radii, Mshells_mass, MfileName, G)
        DfileName = folderName * "/D_t=$t_rounded.txt"
        printToFile(Tshells_radii, Dshells_mass, DfileName, G)
        SfileName = folderName * "/S_t=$t_rounded.txt"
        printToFile(Tshells_radii, Sshells_mass, SfileName, G)
        TfileName = folderName * "/T_t=$t_rounded.txt"
        printToFile(Tshells_radii, Tshells_mass, TfileName, G)
        GPEfileName = folderName * "/GPE_t=$t_rounded"*".txt"
        printToFile_GPE(Tshells_radii, Tshells_GPE, GPEfileName)
    DMTfileName = folderName * "/DMT_t=$t_rounded.txt"
    printToFile(Tshells_radii, Mshells_mass+Dshells_mass, DMTfileName, G)

        timeTaken = (time_ns() - stepStart) / 1e9
        println("Time taken for this step: ", timeTaken, "s\n")
        println(g, t_rounded, "\t", timeTaken)
    end

    # Final results are outputted to a separate folder
    MfileName = folderName_results * "/M_result.txt"
    printToFile(Tshells_radii, Mshells_mass, MfileName, G)
    DfileName = folderName_results * "/D_result.txt"
    printToFile(Tshells_radii, Dshells_mass, DfileName, G)
    SfileName = folderName_results * "/S_result.txt"
    printToFile(Tshells_radii, Sshells_mass, SfileName, G)
    TfileName = folderName_results * "/T_result.txt"
    printToFile(Tshells_radii, Tshells_mass, TfileName, G)
    GPEfileName = folderName_results * "/GPE_result.txt"
    printToFile_GPE(Tshells_radii, Tshells_GPE, GPEfileName)

    timeTaken_total = (time_ns() - functionStart) / 1e9
    println(f, "timeTaken_total=", timeTaken_total)
    println("Total time taken: ", timeTaken_total, "s\n")

    close(f)
    close(g)

    return nothing
end
####################################################################################################

function print_t()
    n=[100,200,1000,2000,3000]

    for i in 1:size(n,1)
        numOfSteps_k = n[i]

        delta_m = (1 - exp(-log(2) * t_end / tau)) / numOfSteps_k  # For having equal amounts of mass decayed every time step
        t_k = zeros(numOfSteps_k + 1)  # Array containing time steps [Gyr]
        t_k[1] = 0
        for j in 2:numOfSteps_k + 1
            t_k[j] = - log(1 - (j - 1) * delta_m) * tau / log(2)
        end

        f=open("time_array_$numOfSteps_k.txt","w")

        for j in 2:size(t_k,1)
            println(f,j,"\t",t_k[j]-t_k[j-1])
        end

        close(f)
    end



end

######################################## Program script ############################################
#print_t()
dmOnly()  # Run the dark matter only algorithm

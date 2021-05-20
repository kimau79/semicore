# using Printf
#using QuadGK

################################## Function that use to solve system of linear equations######################################################

function GPE_extrapolate(r, shells_radii, shells_GPE,a)
#use power law GPE(r) = -Ar^n to extrapolate, for some A, n to be determined
if r>shells_radii[end,3] #the only case that this function should be used
ref_GPE1 = shells_GPE[end-1,1]
ref_GPE2 = shells_GPE[end,1]
ref_radii1 = shells_radii[end-1,3]
ref_radii2 = shells_radii[end,3]
n = log(ref_GPE2/ref_GPE1)/log(ref_radii2/ref_radii1)
A = ref_GPE2/ref_radii2^n
if a==1
return A*r^n	#return normal GPE extrapolation
else
return n*A*r^(n-1)	#return derivative GPE extrapolation
end
end
end


function search_shells_array_index(r,shells_radii)
    
    if r <= shells_radii[1,3]
        return 1
    end

    for i in 1:size(shells_radii,1)-1

        if r > shells_radii[i,3] && r <= shells_radii[i+1,3]

            return i + 1
        end
    end
    return size(shells_radii,1)
end

function gaussian_elimination!(A::Array{Float64,2})

    rows = size(A,1)
    cols = size(A,2)

    # Row index
    row = 1

    # Main loop going through all columns
    for col = 1:(cols-1)

        # finding the maximum element for each column
        max_index = argmax(abs.(A[row:end,col])) + row-1

        # Check to make sure matrix is good!
        if (A[max_index, col] == 0)
            println("matrix is singular!")
            continue
        end

        # swap row with highest value for that column to the top
        temp_vector = A[max_index, :]
        A[max_index, :] = A[row, :]
        A[row, :] = temp_vector

        # Loop for all remaining rows
        for i = (row+1):rows

            # finding fraction
            fraction = A[i,col]/A[row,col]

            # loop through all columns for that row
            for j = (col+1):cols

                 # re-evaluate each element
                 A[i,j] -= A[row,j]*fraction

            end

            # Set lower elements to 0
            A[i,col] = 0
        end
        row += 1
    end
end

function back_substitution(A::Array{Float64,2})

    rows = size(A,1)
    cols = size(A,2)

    # Creating the solution Vector
    soln = zeros(rows)

    for i = rows:-1:1
        sum = 0.0
        for j = rows:-1:i
            sum += soln[j]*A[i,j]
        end
        soln[i] = (A[i, cols] - sum) / A[i, i]
    end

    return soln
end


function gauss_jordan_elimination!(A::Array{Float64,2})

    rows = size(A,1)
    cols = size(A,2)


    # After this, we know what row to start on (r-1)
    # to go back through the matrix
    row = 1
    for col = 1:cols-1
        if (A[row, col] != 0)

            # divide row by pivot and leaving pivot as 1
            for i = cols:-1:col
                A[row,i] /= A[row,col]
            end

            # subtract value from above row and set values above pivot to 0
            for i = 1:row-1
                for j = cols:-1:col
                    A[i,j] -= A[i,col]*A[row,j]
                end
            end
            row += 1
        end
    end
end



######################################### For dmOnly() #############################################
# NFW_params = [rho_0, R_s, c] (see Wiki)
NFW_density(NFW_params, r) = NFW_params[1] / (r / NFW_params[2]) / (1 + r / NFW_params[2]) ^ 2
# NFW_enclosedMass(NFW_params, r) = 4 * pi * NFW_params[1] * NFW_params[2] ^ 3 * (log(1 + r / NFW_params[2]) - r / (NFW_params[2] + r))

# shellRange = [r_1, r_2] where r_1 < r_2
# NFW_shellMass(NFW_params, shellRange) = NFW_enclosedMass(NFW_params, shellRange[2]) - NFW_enclosedMass(NFW_params, shellRange[1])


function NFW_shellMass(NFW_params,shellRange)

    result_integrand(r) = 4 * pi * NFW_params[1] * NFW_params[2] ^ 3 * (NFW_params[2] / (NFW_params[2] + r) + log(NFW_params[2] + r) )

    return result_integrand(shellRange[2]) - result_integrand(shellRange[1])
end

function NFW_EnclosedMass(NFW_params,upperlimit)

    result_integrand(r) = 4 * pi * NFW_params[1] * NFW_params[2] ^ 3 * (NFW_params[2] / (NFW_params[2] + r) + log(NFW_params[2] + r) )

    return result_integrand(upperlimit) - result_integrand(0)
end

function ScaleFactor(Tshells_Vc,Tshells_radii) 
    
    dydx = 100
    
    for i in 2:size(Tshells_radii,1)-1
        dydx = (log(Tshells_Vc[i+1])-log(Tshells_Vc[i-1])) /( (log(Tshells_radii[i+1,2]) - log(Tshells_radii[i-1,2])))
        
        if dydx < 0.3
            Vc = (Tshells_Vc[i-1] + Tshells_Vc[i])/2
            R = (Tshells_radii[i-1,2] + Tshells_radii[i,2])/2
           
            return  Vc, R
        end
        
    end
    println("ScaleFactor error")

end

######################################## For verify_NFW() ##########################################
function NFW_GPE(NFWshells_radii, NFW_params, G)
    NFWshells_GPE = zeros(size(NFWshells_radii, 1))
    for i in 1:size(NFWshells_GPE, 1)
        NFWshells_GPE[i] = -4 * pi * G * NFW_params[1] * NFW_params[2] ^ 3 / NFWshells_radii[i, 3] * log(1 + NFWshells_radii[i, 3] / NFW_params[2])
    end

    return NFWshells_GPE
end


function printToFile_verify_NFW_GPE(fileName, Tshells_radii, Tshells_GPE, NFWshells_GPE,Tshells_enclosedMass,shells_mass)
    f = open(fileName, "w")

    Tshells_Vc = zeros(size(Tshells_radii,1))
    shells_avgRho = zeros(size(Tshells_radii,1))
    shells_rho = zeros(size(Tshells_radii,1))
    for i in 1:size(Tshells_radii, 1)
        Tshells_Vc[i] = sqrt(G * Tshells_enclosedMass[i] / Tshells_radii[i,2]) 
       
        shells_rho[i] = shells_mass[i] / (Tshells_radii[i, 2] ^ 3 - Tshells_radii[i, 1] ^ 3) / (4 / 3 * pi)
        
        shells_avgRho[i] = Tshells_enclosedMass[i] / Tshells_radii[i, 2] ^ 3 / (4 / 3 * pi)
    end
        NFW_Vc_unit, NFW_R_unit = ScaleFactor(Tshells_Vc,Tshells_radii)
    for i in 1:size(Tshells_radii, 1)
        println(f,Tshells_radii[i,2],"\t",Tshells_radii[i,2],"\t",Tshells_radii[i,3],"\t",Tshells_radii[i,2]/NFW_R_unit, "\t", Tshells_Vc[i],"\t", Tshells_Vc[i]/NFW_Vc_unit,"\t",shells_rho[i] ,"\t",shells_avgRho[i],"\t",NFWshells_GPE[i],"\t",Tshells_GPE[i])
    
    end
    close(f)
    return nothing
end



function printToFile_NFW_effPotentialProfile(fileName, Mshells_radii, potentialProfile)
    f = open(fileName, "w")

    for i in 1:size(Mshells_radii, 1)
        println(f, Mshells_radii[i], "\t", potentialProfile[i])
    end
    close(f)
    return nothing
end


####################################################################################################
# Return mass array of NFW profile
# shells_radii = [inner radius, outer radius, shell radius] in the ith row
# shell radius = mid point of boundarys
# shells_mass = [total shell mass]. Assume all mass in a shell concentrate at the position just inside the shell radius
function NFW_shells(NFW_params, numOfShells, shellThicknessFactor,extend_factor)
    #NFW_R_vir = NFW_params[2] * NFW_params[3] 
    
    NFW_R_max = NFW_params[2] * NFW_params[3] * extend_factor

    # Exponentially increasing shellThickness
    #firstShellThickness = NFW_R_vir * (1 - shellThicknessFactor) / (1 - shellThicknessFactor ^ numOfShells)
    firstShellThickness = NFW_R_max * (1 - shellThicknessFactor) / (1 - shellThicknessFactor ^ numOfShells)

    shells_radii = zeros(numOfShells, 3)
    shells_mass = zeros(size(shells_radii, 1))
    for i in 1:size(shells_radii, 1)
        shells_radii[i, 1] = firstShellThickness * (1 - shellThicknessFactor ^ (i - 1)) / (1 - shellThicknessFactor)
        shells_radii[i, 2] = shells_radii[i, 1] + firstShellThickness * shellThicknessFactor ^ (i - 1)
        shells_radii[i, 3] = (shells_radii[i, 1] + shells_radii[i, 2]) / 2
        shells_mass[i] = NFW_shellMass(NFW_params, shells_radii[i, 1:2])
    end

    return shells_radii, shells_mass
end


function enclosedMass(shells_radii, shells_mass)
    #println("EnclosedMass")

    #println("Enclosed mass radius= ",shells_radii[end,1])
    shells_enclosedMass = zeros(size(shells_radii, 1))
    for i in 1:size(shells_enclosedMass, 1)
        shells_enclosedMass[i] = sum(shells_mass[1:i])
    end

    return shells_enclosedMass
end

# Return GPE (per mass) array from a mass array
function GPE(shells_radii, shells_mass, shells_enclosedMass,G)

    
    shells_GPE = zeros(size(shells_radii, 1))
    #=
    shells_rho_appro = zeros(2)
    i_1 = R_lim   
    i_2 = R_lim - 2
    
    shells_rho_appro[1] = shells_mass[i_1] / (shells_radii[i_1, 2] ^ 3 - shells_radii[i_1, 1] ^ 3) / (4 / 3 * pi)
    shells_rho_appro[2] = shells_mass[i_2] / (shells_radii[i_2, 2] ^ 3 - shells_radii[i_2, 1] ^ 3) / (4 / 3 * pi)
    # approximate rho to A r^k
    k = log(shells_rho_appro[1]/shells_rho_appro[2]) / log(shells_radii[i_1,3]/shells_radii[i_2,3])

    A = shells_rho_appro[1] / shells_radii[i_1,3] ^ k
    =#
    
    #println("GPE radius= ",shells_radii[end,1])
    for i in 1:size(shells_radii,1)
        shells_GPE[i] = -G * shells_enclosedMass[i] / shells_radii[i, 3]

        if i < size(shells_radii,1)
            GPEbyOuterShells = 0
            for j in i + 1:size(shells_radii,1)
                GPEbyOuterShells += -G * shells_mass[j] / shells_radii[j, 3]
            end            
            shells_GPE[i] += GPEbyOuterShells
        end
        
    end
    # if power >-2 then we can not obtain a finite analytic GPE value although nfw profile is divergence, the term we add to GPE will be so large
    #=
    if k >=-2
        println("Error on rho approximation: k")
        println("k=",k)
    else
        for i in 1:size(shells_radii,1)
            # GPE by outer shells R_lim to infinity
            shells_GPE[i] +=  4 / (k + 2) * pi * G *  A * shells_radii[end ,3] ^ (k + 2)
        end
    end
    =#
    
    return shells_GPE

end


# Return angular momentum (per mass) array
function L(shells_radii, shells_enclosedMass, G)
    shells_L = zeros(size(shells_radii, 1))
    for i in 1:size(shells_L, 1)
        shells_L[i] = (G * shells_enclosedMass[i] * shells_radii[i, 3]) ^ (1 / 2)
    end

    return shells_L
end


# Return total energy (per mass) array of any just-decayed particle at different radii
function totalE_afterDecay(shells_radii, shells_GPE, shells_L, v_k)
    shells_totalE_afterDecay = zeros(size(shells_radii, 1))
    for i in 1:size(shells_totalE_afterDecay, 1)
        shells_totalE_afterDecay[i] = shells_GPE[i] + (shells_L[i] / shells_radii[i, 3]) ^ 2 / 2 + v_k ^ 2 / 2
    end

    return shells_totalE_afterDecay
end


function energyEquation(r, L, totalE_afterDecay, Tshells_radii, Tshells_GPE)
    if r <= 0  # Rejected

        return zeros(NaN)  # Error
    elseif r <= Tshells_radii[1, 3]  # r small
        return Tshells_GPE[1] + (L / r) ^ 2 / 2 - totalE_afterDecay
    elseif r > Tshells_radii[end,3]
        return  GPE_extrapolate(r,Tshells_radii, Tshells_GPE,1)  + (L / r) ^ 2 / 2 - totalE_afterDecay  
    else  # r in between; value by interpolation
        radiusIndex = -1  # Just for the definition
        for i in 2:size(Tshells_radii, 1)
            if r <= Tshells_radii[i, 3]
                radiusIndex = i
            
                break
            end
        end
        intervalSlope = (Tshells_GPE[radiusIndex] - Tshells_GPE[radiusIndex - 1]) / (Tshells_radii[radiusIndex, 3] - Tshells_radii[radiusIndex - 1, 3])
        intervalIntercept = Tshells_GPE[radiusIndex] - intervalSlope * Tshells_radii[radiusIndex, 3]
        radiusGPE = intervalSlope * r + intervalIntercept

        return radiusGPE + (L / r) ^ 2 / 2 - totalE_afterDecay
    end
end


function ellipseSolver(r_0, L, totalE_afterDecay, shells_radii, Tshells_GPE, tol_ellipseGuess)
    # Search in [l1, l2] U [r1, r2] using the bisection method

    firstShellThickness = shells_radii[1, 2]  # To be used as a tolerance
   
    # Some initial checking
    if energyEquation(r_0, L, totalE_afterDecay, shells_radii, Tshells_GPE) >= 0
        # This should not happen unless GPE/totalE are not updated properly (= 0 occurs when v_k = 0)
        println("ellipseSolver: v_k probably too small; no solvable roots")
        
        # println(energyEquation(r_0, L, totalE_afterDecay, Tshells_radii, Tshells_GPE, Tshells_enclosedMass))
        # zeros(NaN)  # Halt program
        return r_0, r_0  # If this happens, radii just stay put (i.e. solution for v_k = 0)
    elseif totalE_afterDecay >= 0  # Escaped
        return -1, -1
    else  # If checking passed
        l2 = r_0
        r1 = r_0
    end
    #println("start r_0 =",r_0)
    # Setting l1 and r2
    l1 = firstShellThickness
    while energyEquation(l1, L, totalE_afterDecay, shells_radii, Tshells_GPE) <= 0
        l1 /= 2
    end
    #println("l1 end")
    r2 = shells_radii[end, 3]
    while energyEquation(r2, L, totalE_afterDecay, shells_radii, Tshells_GPE) <= 0
        r2 *= 2
    end
    #println("r2 end")
    # Bisection method
    lastDiff = 0
    while (l2 - l1 > firstShellThickness * tol_ellipseGuess*0.001) && (l2 - l1 != lastDiff)
        lastDiff = l2 - l1
        l3 = (l1 + l2) / 2
        energyEquation_value = energyEquation(l3, L, totalE_afterDecay, shells_radii, Tshells_GPE)
        if energyEquation_value < 0
            l2 = l3
        elseif energyEquation_value > 0
            l1 = l3
        else
            l1 = l3
            l2 = l3
        end
    end
    lastDiff = 0
    while (r2 - r1 > firstShellThickness * tol_ellipseGuess*0.001) && (r2 - r1 != lastDiff)
        lastDiff = r2 - r1
        r3 = (r2 + r1) / 2
        energyEquation_value = energyEquation(r3, L, totalE_afterDecay, shells_radii, Tshells_GPE)
        if energyEquation_value < 0
            r1 = r3
        elseif energyEquation_value > 0
            r2 = r3
        else
            r1 = r3
            r2 = r3
        end
    end

    root1 = (l1 + l2) / 2
    root2 = (r1 + r2) / 2
    return root1, root2
end

# Return ellipse array for two roots

function ellipseRadii(shells_L, shells_totalE_afterDecay, Tshells_radii, Tshells_GPE, tol_ellipseGuess)
    shells_ellipseRadii = zeros(size(Tshells_radii, 1), 2)

    for i in 1:size(shells_ellipseRadii, 1)
        
        root1, root2 = ellipseSolver(Tshells_radii[i, 3], shells_L[i], shells_totalE_afterDecay[i], Tshells_radii, Tshells_GPE, tol_ellipseGuess)

        shells_ellipseRadii[i, 1] = root1
        shells_ellipseRadii[i, 2] = root2
    end

    return shells_ellipseRadii
end


function weightFactorArray(r_ref, shells_ellipseRadii, L, shells_totalE, Tshells_GPE, Tshell_radii, Tshells_enclosedMass, t_i, orderOfpolynomial, G, NFW_params)
    weightFactor = zeros(size(shells_ellipseRadii, 1))
   

    for i in 1:size(weightFactor, 1)  # Looping each r_0
        r_max = shells_ellipseRadii[i, 2]
        r_min = shells_ellipseRadii[i, 1]

        if r_max == -1 && r_min == -1  # Escaped the whole system
         
            weightFactor[i] = 0
        elseif r_min > r_ref
            weightFactor[i] = 0
        elseif r_max <= r_ref
            weightFactor[i] = 1
        else
            weightFactor[i] = weightFactorSolver(r_ref, r_max, r_min, L[i], shells_totalE[i], Tshells_GPE, Tshell_radii, Tshells_enclosedMass, t_i, orderOfpolynomial, G, NFW_params)

        end
    end
    return weightFactor
end

function U_eff_gfunction(x,r_max,r_min, L, Tshells_radii, Tshells_GPE)
    # convert r to be x between (0,1)
    r = x * (r_max - r_min) +r_min

    if r <= 0  # Rejected

        return zeros(NaN)  # Error
    elseif r <= Tshells_radii[1, 3]  # r small
        return Tshells_GPE[1]  + (L / r) ^ 2 / 2 
    elseif r > Tshells_radii[end, 3]  # r big
        
        return GPE_extrapolate(r,Tshells_radii, Tshells_GPE,1) + (L / r) ^ 2 / 2
    else  # r in between; value by interpolation
        radiusIndex = -1  # Just for the definition
        for i in 2:size(Tshells_radii, 1)
            if r <= Tshells_radii[i, 3]
                radiusIndex = i
                break
            end
        end
        intervalSlope = (Tshells_GPE[radiusIndex] - Tshells_GPE[radiusIndex - 1]) / (Tshells_radii[radiusIndex, 3] - Tshells_radii[radiusIndex - 1, 3])
        intervalIntercept = Tshells_GPE[radiusIndex] - intervalSlope * Tshells_radii[radiusIndex, 3]
        radiusGPE = intervalSlope * r + intervalIntercept

        return radiusGPE + (L / r) ^ 2 / 2 
    end
end

function dU_eff(r, L, Tshells_radii, Tshells_enclosedMass, Tshells_GPE)
    if r <= 0  # Rejected
        println("Dead End")
        return zeros(NaN)  # Error
    elseif r <= Tshells_radii[1, 3]  # r small
        return G * Tshells_enclosedMass[1] / Tshells_radii[1,3] / r   - (L / r) ^ 2 / r
    elseif r > Tshells_radii[end, 3]  # r big
        return GPE_extrapolate(r,Tshells_radii, Tshells_GPE,2)  - (L / r) ^ 2 / r
    else  # r in between; value by interpolation
        radiusIndex = -1  # Just for the definition
        for i in 2:size(Tshells_radii, 1)
            if r <= Tshells_radii[i, 3]
                radiusIndex = i
                break
            end
        end
        intervalSlope = (Tshells_enclosedMass[radiusIndex] - Tshells_enclosedMass[radiusIndex - 1]) / (Tshells_radii[radiusIndex, 3] - Tshells_radii[radiusIndex - 1, 3])
        intervalIntercept = Tshells_enclosedMass[radiusIndex] - intervalSlope * Tshells_radii[radiusIndex, 3]
        radiusMass = intervalSlope * r + intervalIntercept

        return G * radiusMass / r / r  - (L / r) ^ 2 / r
    end

end

function dU_eff_NFW(r, NFW_params, G,L)
    R_s = NFW_params[2]  
    rho_0 = NFW_params[1]

    dUdr = 4 * pi * G * rho_0 * R_s ^ 3 / r ^ 2 * log(1 + r / R_s) - 4 * pi * G * rho_0 * R_s ^ 2 / r / (1 + r / R_s)   

    return dUdr - L ^ 2 /  r ^ 3

end
###### Approximate using straight line y=mx+c for [0,0.02] 
function res_StraightLine(res, res_node_0)

    # x_node for compute straight line approximation on res(x)
    x_node = [0,0.002,0.01,0.02]
    ## Compute the value of res(x) at each x_node
    res_node = zeros(size(x_node,1))
    ### slope = m , constant_c = c
    slope = zeros(size(x_node,1) - 1)
    constant_c = zeros(size(x_node,1) - 1)
    ### I_res[i] reprsent the result of integrate res(x) from x_node[i] to x_node[i+1]
    Ires_sl= zeros(size(x_node,1) - 1)
    #### compute the first node of res which has strong singularity if we calculate it computationally
    res_node[1] = res_node_0
    ##### compute the reset of the nodes
    for i in 2:size(x_node,1)
        res_node[i] = res(x_node[i])
    end
    ### compute slope = m , constant_c = c
    for i in 1:size(slope,1)
        slope[i] = (res_node[i + 1] - res_node[i]) /(x_node[i + 1] - x_node[i])
        constant_c[i] = res_node[i] - slope[i] * x_node[i]
    end
    Ires_sl_max = 0

    for i in 1:size(Ires_sl,1)
        f(x) =  1/2 * slope[i] * x ^ 2 + constant_c[i] * x 

        Ires_sl[i] = f(x_node[i + 1]) - f(x_node[i]) 
        Ires_sl_max += Ires_sl[i]
    end

    return x_node, slope, constant_c, Ires_sl, Ires_sl_max

end
###### Approximate using polynomial y= a +  b x + c x^2 + d x^3 ... for [0.02,1] 
function res_Polynomial(res, x_a, Number_of_point, res_node_1, x2, x1)
    
    # x_node for compute straight line approximation on res(x)
    x_node = zeros(Number_of_point)
    ## x_0 = first x node for polynomial approximation
    dx = (x2 - x1) / (Number_of_point - 1)
    for i in 1:Number_of_point
        x_node[i] = x1 + (i - 1) * dx
    end
    
    ## Compute the value of res(x) at each x_node
    res_node = zeros(Number_of_point)
    for i in 1:Number_of_point - 1
       res_node[i] = res(x_node[i])
    end
    res_node[end] = res_node_1

    # set up linear euqation for approximating res(x) numerically

    # res(x) = a +  b x + c x^2 + d x^3 ...up to x^n order ,n = number of points -1 
     
    equation_matrix = zeros(Number_of_point,Number_of_point + 1)

    for i in 1:Number_of_point 
        # choose x_node[i] for i-th linear equation 
        for j in 1:Number_of_point
            equation_matrix[i,j] = (x_node[i] - x_a) ^ (j - 1)
        end
        equation_matrix[i,end] = res_node[i]
    end
    
    ## solution for the linear equation 
    soln = zeros(Number_of_point)
    
    # solve the matrix by gaussian elimination 
    gaussian_elimination!(equation_matrix)
    gauss_jordan_elimination!(equation_matrix)
    soln = back_substitution(equation_matrix)   

    return soln    


end

function weightFactorSolver(r_ref,r_max,r_min,L,E,Tshells_GPE,Tshells_radii,Tshells_enclosedMass,t_i,orderOfpolynomial,G,NFW_params)

    U_eff(x) = U_eff_gfunction(x,r_max,r_min,L,Tshells_radii,Tshells_GPE)
  
    # convert the variable r to x between (0,1) so that the it can be normalized
    x_min = 0  
    x_max = 1
    x_ref = (r_ref - r_min) / (r_max - r_min)
    
    if t_i == 2
        dU_max = dU_eff_NFW(r_max,NFW_params,G,L) *(r_max - r_min)
        dU_min = dU_eff_NFW(r_min,NFW_params,G,L) *(r_max - r_min)
    else
        dU_max = dU_eff(r_max,L,Tshells_radii,Tshells_enclosedMass,Tshells_GPE) * (r_max - r_min)
        dU_min = dU_eff(r_min,L,Tshells_radii,Tshells_enclosedMass,Tshells_GPE) * (r_max - r_min)
    end

    res(x)= 1 / sqrt( E - U_eff(x) )- 1 / sqrt( -1 * dU_min * (x - x_min) )- 1 / sqrt( dU_max * (x_max - x) )
     #### compute the first and final node of res which has strong singularity if we calculate it computationally
     res_node_0 = -1 / sqrt(dU_max * (x_max - x_min))
     res_node_1 = -1 / sqrt(-dU_min * (x_max -x_min))
       
    ####### Function for computing size(sl_xnode,1) - 1 straight line for res(x) at [0,0.2]
    x_node, slope, constant_c, Ires_sl, Ires_sl_max = res_StraightLine(res, res_node_0)
    # reference point for taylor series expansion (res(x) approximation)
    x_a = (x_max - x_node[end]) / 2
    Number_of_point = orderOfpolynomial + 1
    ####### Function for computing x-order Polynomial for res at [0.2,1]
    soln = res_Polynomial(res, x_a, Number_of_point, res_node_1, x_max, x_node[end])
    ###########################################
    Ires_poly_max = Ires_poly_min = 0
    ### Compute I_res_max = integrate res(x)dx from 0 to 1
    for i in 1:Number_of_point

        Ires_poly_max += soln[i] * (x_max - x_a) ^ (i) /i
        Ires_poly_min += soln[i] * (x_node[end] -  x_a) ^ (i) /i

    end 

    Ires_max = Ires_sl_max + Ires_poly_max - Ires_poly_min
    
    ######################################## r_ref in [0,0.002]#######################################
    Ires_ref = 0 

    # Compute I_res_ref = integrate res(x)dx from 0 to x_ref
    if x_ref < x_node[end]
        x_indicator = 1 
        while x_ref > x_node[x_indicator + 1]
    
            if x_indicator + 1 == size(x_node,1)
                
                break
            end 
            x_indicator += 1
                
        end
        #####
        Ires_f(x) = 1 / 2 * slope[x_indicator] * x ^ 2 + constant_c[x_indicator] * x 
        # app y = mx + c 
    
        for i in 1:x_indicator - 1 
            Ires_ref += Ires_sl[i]
        end
        Ires_ref += Ires_f(x_ref) - Ires_f(x_node[x_indicator])
    else 
        Ires_ref += Ires_sl_max
        ########################################## r_ref in [0.002,1]##########################################################
        for i in 1:Number_of_point
            Ires_ref += soln[i] *(x_ref - x_a) ^ (i) / i
        end
        Ires_ref -= Ires_poly_min 
    end
    
    ##################################################
    # 1 / sqrt(E- U_eff(x)) = res(x) + res(x) + 1 / sqrt( -1 * dU_min * (x - x_min) ) + 1 / sqrt( dU_max * (x_max - x) )
    # N1 and N2 are analytic solution of integrae 1 / sqrt( -1 * dU_min * (x - x_min) ) and 1 / sqrt( dU_max * (x_max - x) ) from x_min to x_ref 
    # D1 and D2 are analytic solution of integrae 1 / sqrt( -1 * dU_min * (x - x_min) ) and 1 / sqrt( dU_max * (x_max - x) ) from x_min to x_max 
    
    N1 = 2 * sqrt( (x_ref - x_min) / (-dU_min ) )
    N2 = 2 * ( sqrt(x_max - x_min) - sqrt(x_max - x_ref) ) / sqrt(dU_max)

    D1 = 2 * sqrt( (x_max - x_min) / (-dU_min ) )
    D2 = 2 * sqrt( (x_max - x_min) / (dU_max ) )

    nominator = Ires_ref + N1 + N2 
    denominator = Ires_max + D1 + D2


    return nominator/denominator
end


function updateShellsMass(shells_radii, shells_ellipseRadii, Mshells_mass, p_undecayed,L,shells_totalE,Tshells_GPE,Tshells_enclosedMass,t_i,orderOfpolynomial, G, NFW_params)
    Mshells_decayedMass = Mshells_mass * (1 - p_undecayed)  # To be redistributed
    Mshells_mass *= p_undecayed  # Remaining mass

    Dshells_enclosedMass_decayedMass = zeros(size(shells_radii, 1))
    
    for i in 1:size(Dshells_enclosedMass_decayedMass, 1)

        weightFactor = weightFactorArray(shells_radii[i, 2],shells_ellipseRadii,L,shells_totalE,Tshells_GPE,shells_radii,Tshells_enclosedMass,t_i, orderOfpolynomial, G, NFW_params)
        Dshells_enclosedMass_decayedMass[i] = sum(Mshells_decayedMass .* weightFactor)
        
    end

    Dshells_decayedMass = zeros(size(Dshells_enclosedMass_decayedMass, 1))
    if Dshells_decayedMass != []  # If all mothers at all radius escape upon decay
        Dshells_decayedMass[1] = Dshells_enclosedMass_decayedMass[1]
        for i in 2:size(Dshells_decayedMass, 1)
            Dshells_decayedMass[i] = Dshells_enclosedMass_decayedMass[i] - Dshells_enclosedMass_decayedMass[i - 1]
        end
    end



    return Mshells_mass, Dshells_decayedMass
end

function adiabaticExpansion(shells_radii, shells_mass, Tshells_enclosedMass, Tshells_enclosedMass_updated)
    # At this moment:
    # Mshells_radii is short as original
    # Dshells_radii is extended
    # Tshells_radii is short as original
    # Tshells_radii_updated is extended

    # if size(Tshells_enclosedMass, 1) < size(shells_radii, 1)
    #     println("adiabaticExpansion: shell sizes do not match by", size(shells_radii, 1) - size(Tshells_enclosedMass, 1))
    #     for i in 1:size(shells_radii, 1) - size(Tshells_enclosedMass, 1)
    #         push!(Tshells_enclosedMass, Tshells_enclosedMass[end])
    #     end
    # end
    
    expansionRatios = Tshells_enclosedMass[1:size(shells_radii, 1)] ./ Tshells_enclosedMass_updated[1:size(shells_radii, 1)]

ER_filename = folderName
run_id2 = -1
for i=1:100
if !isfile(ER_filename*"/ER"*string(i)*".txt")
	run_id2 = i
	break
end
end
ER_filename *="/ER"*string(run_id2)*".txt"
writedlm(ER_filename,expansionRatios)
 
    # # Hotfix for expansion ratio very close to 1 (maybe not)
    # for i in 1:size(expansionRatios, 1)
    #     expansionRatios[i] = round(expansionRatios[i], digits=3)  # I just picked digits=3
    # end
    # To check if it actaully contracts instead of expanding. But this doesn't really matter
    contractionCount = count(i -> (i < 1), expansionRatios)
    if contractionCount > 0
        # println("adiabaticExpansion: expansion ratio smaller than 1, i.e. NOT expanding. Count: ", contractionCount, ", min ratio: ", findmin(expansionRatios)[1])
        # zeros(NaN)  # To cause error, halting the program
    end

    # shells_expandedRadii = shells_radii[:, 3] .* expansionRatios
    shells_expandedRadii = shells_radii[:, 2] .* expansionRatios  # Use 2 or 3? 2

    # To make sure expandedRadii is "monotonic" (never seen useful)
    violationCount = 0
    checkedEntry = -1
    while checkedEntry != size(shells_expandedRadii, 1) - 1
        checkedEntry = -1
        for i in 1:size(shells_expandedRadii, 1) - 1
            if shells_expandedRadii[i] > shells_expandedRadii[i + 1]
                violationCount += 1

                eR_1 = shells_expandedRadii[i]
                eR_2 = shells_expandedRadii[i + 1]
                shells_expandedRadii[i] = eR_2
                shells_expandedRadii[i + 1] = eR_1

                break
            else
                checkedEntry = i
            end
        end
    end
    
    if violationCount > 0
        println("adiabaticExpansion: violationCount = ", violationCount)
    end

    
    expandedShells_radii = shells_radii
    
    expandedShells_mass = zeros(size(shells_radii,1),1)
    
    for i in 1:size(shells_radii, 1)  # This interpolation thing should work if the relation is monotonic. Check total mass after expansion.
        e1 = shells_radii[i, 1]  # Inner radius of expanded shells
        e2 = shells_radii[i, 2]  # Outer radius of expanded shells
        
        e1_smallerThanID = -1
        for j in 1:size(shells_expandedRadii, 1)
            if e1 < shells_expandedRadii[j]
                e1_smallerThanID = j
                break
            end
        end
        
        e2_smallerThanID = -1
        for j in 1:size(shells_expandedRadii, 1)
            if e2 < shells_expandedRadii[j]
                e2_smallerThanID = j
                break
            end
        end
        
        
        
        if e1_smallerThanID == 1
            m = (shells_radii[e1_smallerThanID, 2] - 0) / (shells_expandedRadii[e1_smallerThanID] - 0)
            c = 0
            r1 = m * e1 + c
        elseif e1_smallerThanID != -1
            m = (shells_radii[e1_smallerThanID, 2] - shells_radii[e1_smallerThanID - 1, 2]) / (shells_expandedRadii[e1_smallerThanID] - shells_expandedRadii[e1_smallerThanID - 1])
            c = shells_radii[e1_smallerThanID, 2] - m * shells_expandedRadii[e1_smallerThanID]
            r1 = m * e1 + c
        else
                        
            r1 = -1  # Should never happen
        end
        
        if e2_smallerThanID == 1
            m = (shells_radii[e2_smallerThanID, 2] - 0) / (shells_expandedRadii[e2_smallerThanID] - 0)
            c = 0
            r2 = m * e2 + c
        elseif e2_smallerThanID != -1
            m = (shells_radii[e2_smallerThanID, 2] - shells_radii[e2_smallerThanID - 1, 2]) / (shells_expandedRadii[e2_smallerThanID] - shells_expandedRadii[e2_smallerThanID - 1])
            c = shells_radii[e2_smallerThanID, 2] - m * shells_expandedRadii[e2_smallerThanID]
            r2 = m * e2 + c
        else
            r2 = -1  # Will happen once
            # println("adiabaticExpansion: r2 = -1")
        end
        
        firstShellThickness = shells_radii[1, 2]
        shellThicknessFactor = (shells_radii[2, 2] - shells_radii[2, 1]) / firstShellThickness
        if r1 != -1
            totalLen = 0
            r1_smallerThanID = 0
            while totalLen <= r1
                r1_smallerThanID += 1
                totalLen += firstShellThickness * shellThicknessFactor ^ (r1_smallerThanID - 1)
            end
            
            if r1_smallerThanID > size(shells_radii, 1)
                println("adiabatic Expansion error: r1 > outermost radius")
                continue  # Hotfix to weird boundary cases
            end
            
        else
            
            println("adiabaticExpansion error: r1 = -1")  # Prompt error
           
            continue  # Hotfix to weird boundary cases
        end
        
        if r2 != -1
            totalLen = 0
            r2_smallerThanID = 0
            while totalLen <= r2
                r2_smallerThanID += 1
                totalLen += firstShellThickness * shellThicknessFactor ^ (r2_smallerThanID - 1)
            end
        else
            r2_smallerThanID = -1  # Special treatment
        end
       
       
        expandedShells_mass[i] += shells_mass[r1_smallerThanID] * (1 - (r1 ^ 3 - shells_radii[r1_smallerThanID, 1] ^ 3) / (shells_radii[r1_smallerThanID, 2] ^ 3 - shells_radii[r1_smallerThanID, 1] ^ 3))
        if r2_smallerThanID == -1
            expandedShells_mass[i] += shells_mass[end]  # This is why the density is always weird at the end (solved)
            r2_smallerThanID = size(shells_radii, 1)
        else
            
            expandedShells_mass[i] += shells_mass[r2_smallerThanID] * (1 - (shells_radii[r2_smallerThanID, 2] ^ 3 - r2 ^ 3) / (shells_radii[r2_smallerThanID, 2] ^ 3 - shells_radii[r2_smallerThanID, 1] ^ 3))
        end
        
        if r1_smallerThanID == r2_smallerThanID
            expandedShells_mass[i] -= shells_mass[r1_smallerThanID]
        elseif r2_smallerThanID - r1_smallerThanID > 1
            expandedShells_mass[i] += sum(shells_mass[r1_smallerThanID + 1:r2_smallerThanID - 1])
        end
        
    end

    return expandedShells_mass
end


function printToFile(shells_radii, shells_mass, fileName, G)
    f = open(fileName, "w")
    
    shells_rho = zeros(size(shells_radii, 1))  # Shell density
    shells_enclosedMass = zeros(size(shells_radii, 1))  # Enclosed mass
    shells_avgRho = zeros(size(shells_radii, 1))  # Average density
    shells_Vcir = zeros(size(shells_radii, 1))  # Circular velocity
    for i in 1:size(shells_rho, 1)
        shells_rho[i] = shells_mass[i] / (shells_radii[i, 2] ^ 3 - shells_radii[i, 1] ^ 3) / (4 / 3 * pi)
        shells_enclosedMass[i] = sum(shells_mass[1:i])
        shells_avgRho[i] = shells_enclosedMass[i] / shells_radii[i, 2] ^ 3 / (4 / 3 * pi)
        shells_Vcir[i] = sqrt(G * shells_enclosedMass[i] / shells_radii[i, 2]) 
    end

    for i in 1:size(shells_radii, 1)
        println(f, shells_radii[i, 1], "\t", shells_radii[i, 2], "\t", shells_radii[i, 3], "\t", shells_mass[i], "\t", shells_rho[i], "\t", shells_enclosedMass[i], "\t", shells_avgRho[i], "\t", shells_Vcir[i])
    end
    
    close(f)
    return nothing
end


function printToFile_GPE(Tshells_radii, Tshells_GPE, fileName)
    f = open(fileName, "w")

   
    
    for i in 1:size(Tshells_radii, 1)
        println(f, Tshells_radii[i, 1], "\t", Tshells_radii[i, 2], "\t", Tshells_radii[i, 3], "\t", Tshells_GPE[i])
    end
    close(f)
    return nothing
end




######################################### Not used yet #############################################
function shellTrimmer(shells_radii, shells_mass)
    numOfZeros = 0
    for i in 0:size(shells_radii, 1) - 1
        if shells_mass[end - i] == 0
            numOfZeros += 1
        else
            break
        end
    end

    return shells_radii[1:end - numOfZeros, :], shells_mass[1:end - numOfZeros]
end


# Removing the "Boltzmann tail" of baryon particles (error: T should be r dependent)
function barEscape(T, Tshells_GPE, Bshells_mass, m, k)
    totalBarMass = sum(Bshells_mass)

    Bshells_escapeV = zeros(size(Bshells_mass, 1))
    for i in 1:size(Bshells_escapeV, 1)
        Bshells_escapeV[i] = (-Tshells_GPE[i] * 2) ^ (1 / 2)
    end

    integrand(v) = 4 * pi * v ^ 2 * (m / (2 * pi * k * T)) ^ (3 / 2) * exp(-m * v ^ 2 / (2 * k * T))
    for i in 1:size(Bshells_mass, 1)
        retainedFraction = quadgk(integrand, 0, Bshells_escapeV[i])[1]
        Bshells_mass[i] *= retainedFraction
    end

    totalBarMass_updated = sum(Bshells_mass)

    # println("v_rms = ", (3 * k * T / m) ^ (1 / 2), " kpc / s, max escapeV = ", findmax(Bshells_escapeV)[1], " kpc / s, min escapeV = ", findmin(Bshells_escapeV)[1], " kpc / s")
    println("% escaped: ", (1 - totalBarMass_updated / totalBarMass) * 100, "%")

    return totalBarMass_updated
end


function escapedRemoval(Tshells_enclosedMass, Tshells_GPE_updated, shells_radii, shells_mass, G)
    for i in 1:size(shells_radii, 1)
        KE = G * Tshells_enclosedMass[i] / (2 * shells_radii[i, 3])  # Assume circularly moving particles
        # println(Tshells_GPE_updated[i], "\t", KE)
        if Tshells_GPE_updated[i] + KE > 0
            shells_mass[i] = 0
        end
    end

    shells_radii, shells_mass = shellTrimmer(shells_radii, shells_mass)
    return shells_mass
end


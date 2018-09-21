if !isdefined(:ExoplanetsSysSim) using ExoplanetsSysSim end
using Distributions
using Roots #if we want to use the 'find_zero()' function
using Optim

##### This code is a translation of 'MRpredict.R' at 'https://github.com/Bo-Ning/Predicting-exoplanet-mass-and-radius-relationship/blob/master/MR-predict/MRpredict.R'

struct MR_param_Ning2018
    Mass_min::Float64 #log10 scale
    Mass_max::Float64 #log10 scale
    Radius_min::Float64 #log10 scale
    Radius_max::Float64 #log10 scale
    deg::Int64
    weights_mle::Vector{Float64}
    indices_keep::Vector{Int64}
    indices_div::Vector{Int64}
    indices_rem::Vector{Int64}
end

function pdfnorm_beta(x::Float64, x_obs::Float64, x_sd::Float64, x_max::Float64, x_min::Float64, shape1::Int64, shape2::Int64 ; log::Bool=true)
    if log == true
        norm_dist = Normal(10^x, x_sd)
    else
        norm_dist = Normal(x, x_sd)
    end
    #beta_dist = Beta(shape1, shape2)
    return pdf.(norm_dist, x_obs) .* pdf.(Beta(shape1, shape2), (x-x_min)/(x_max-x_min)) / (x_max-x_min)
end

function fn_for_integrate(x_obs::Float64, x_sd::Float64, deg::Int, degree::Int, x_max::Float64, x_min::Float64 ; log::Bool=false, abs_tol::Float64=1e-10)
    shape1, shape2 = degree, deg-degree+1
    return quadgk(x -> pdfnorm_beta(x, x_obs, x_sd, x_max, x_min, shape1, shape2; log=log), x_min, x_max, abstol=abs_tol)[1]
end

function conditional_density(y::Float64, y_max::Float64, y_min::Float64, x_max::Float64, x_min::Float64, deg::Int64, w_hat::Vector{Float64}, indices_keep::Vector{Int64}, indices_div::Vector{Int64}, indices_rem::Vector{Int64}; y_sd::Union{Void,Float64}=nothing, qtl::Vector{Float64}=[0.16, 0.84])
    deg_vec = 1:deg

    #To compute conditional mean, variance, quantile, distribution:
    y_beta_indv = zeros(length(indices_rem))
    if y_sd == nothing
        #y_stdardize = (y-y_min)/(y_max-y_min)
        #y_beta_indv = map(degree -> pdf.(Beta(degree, deg-degree+1), (y-y_min)/(y_max-y_min)), deg_vec)::Array{Float64,1} / (y_max-y_min)
        map!(degree -> pdf.(Beta(degree, deg-degree+1), (y-y_min)/(y_max-y_min)), y_beta_indv, deg_vec[indices_rem])
        y_beta_indv /= (y_max-y_min)
    else
        #y_beta_indv = map(degree -> fn_for_integrate(y, y_sd, deg, degree, y_max, y_min), deg_vec)::Array{Float64,1}
        map!(degree -> fn_for_integrate(y, y_sd, deg, degree, y_max, y_min), y_beta_indv, deg_vec[indices_rem])
    end
    #y_beta_pdf = kron(ones(deg), y_beta_indv)[indices_keep]
    #y_beta_pdf = y_beta_indv[indices_rem]
    denominator = sum(w_hat .* y_beta_indv)

    #y_beta_indv = y_beta_indv[indices_rem]

    #Mean:
#mean_beta_indv = deg_vec/(deg+1)*(x_max-x_min)+x_min
    #mean_beta = kron(mean_beta_indv, y_beta_indv)[indices_keep]
#mean_beta = mean_beta_indv[indices_div] .* y_beta_indv
#mean_numerator = sum(w_hat .* mean_beta::Vector{Float64})
#mean = mean_numerator / denominator
    mean = sum(w_hat::Vector{Float64} .* (((deg_vec::UnitRange{Int64})/(deg+1)*(x_max-x_min)+x_min)[indices_div] .* y_beta_indv::Vector{Float64})::Vector{Float64}) / denominator::Float64
    #mean = sum(w_hat .* kron(deg_vec/(deg+1)*(x_max-x_min)+x_min, y_beta_indv)) / denominator

    #Variance:
#var_beta_indv = deg_vec .* (deg-deg_vec+1) / ((deg+1)^2*(deg+2))*(x_max-x_min)^2
    #var_beta = kron(var_beta_indv, y_beta_indv)[indices_keep]
#var_beta = var_beta_indv[indices_div] .* y_beta_indv
#var_numerator = sum(w_hat .* var_beta::Vector{Float64})
#var = var_numerator /denominator
    var = sum(w_hat::Vector{Float64} .* (((deg_vec::UnitRange{Int64}) .* (deg-deg_vec+1) / ((deg+1)^2*(deg+2))*(x_max-x_min)^2)[indices_div] .* y_beta_indv::Vector{Float64})::Vector{Float64}) / denominator::Float64
    #var = sum(w_hat .* kron(deg_vec .* (deg-deg_vec+1) / ((deg+1)^2*(deg+2))*(x_max-x_min)^2, y_beta_indv)) / denominator

    #Quantile:
    function pbeta_conditional_density(x::Float64)::Float64
        function mix_density(j::Float64)::Float64
            #x_indv_cdf::Vector{Float64} = map(degree -> cdf.(Beta(degree, deg-degree+1), (j-x_min)/(x_max-x_min)), deg_vec)
            ###quantile_beta = x_indv_cdf[indices_div] .* y_beta_indv[indices_rem]
            #quantile_numerator = sum(w_hat .* kron(x_indv_cdf, y_beta_indv)::Vector{Float64})
            #return quantile_numerator/denominator
            #return sum(w_hat .* kron(map(degree -> cdf.(Beta(degree, deg-degree+1), (j-x_min)/(x_max-x_min)), deg_vec), y_beta_indv::Vector{Float64})[indices_keep]) / denominator
            @fastmath @inbounds @views return sum(w_hat::Vector{Float64} .* (map(degree -> cdf.(Beta(degree, deg::Int64-degree+1), (j-x_min)/(x_max-x_min)), (deg_vec::UnitRange{Int64})[indices_div]) .* y_beta_indv::Vector{Float64})) / denominator::Float64
        end

        return mix_density(x)
    end

    function conditional_quantile(q::Float64)::Float64
        function g(x::Float64)::Float64
            return pbeta_conditional_density(x) - q
        end
        function g2(x::Float64)::Float64
            return (pbeta_conditional_density(x) - q)^2
        end

        #return find_zero(g, (x_min, x_max), FalsePosition())
        #return find_zero(g, (x_min, x_max))::Float64 ###'find_zero' is from 'using Roots' (is this a good function to use?); Eric suggested brent's method in the Optim package

        result = optimize(g2, x_min, x_max, Brent(), rel_tol=1e-3, abs_tol=1e-3) #Brent(); GoldenSection()
        if Optim.converged(result)
            return Optim.minimizer(result)::Float64
        else
            return find_zero(g, (x_min, x_max))::Float64
        end
    end

    quantile = map(conditional_quantile, qtl)

    return (mean, var, quantile, denominator, y_beta_indv)
end

function predict_mass_given_radius(Radius::Float64, param::MR_param_Ning2018 ; R_sigma::Union{Void,Float64}=nothing, posterior_sample::Bool=false, qtl::Vector{Float64}=[0.16, 0.84])
    #Convert data to log scale:
    l_radius = log10.(Radius)

    #The 'posterior_sample == false' condition can deal with two cases:
    #Case I: if input data do not have measurement errors
    #Case II: if input data have measurement error
    if posterior_sample == false
        predicted_value = conditional_density(l_radius, param.Radius_max, param.Radius_min, param.Mass_max, param.Mass_min, param.deg, param.weights_mle, param.indices_keep, param.indices_div, param.indices_rem; y_sd=R_sigma, qtl=qtl)
        predicted_mean = predicted_value[1]
        predicted_quantiles = predicted_value[3]
        #predicted_lower_quantile = predicted_value[3][1]
        #predicted_upper_quantile = predicted_value[3][2]
    elseif posterior_sample == true
        #Case III: if the input are posterior samples
        radius_sample = log10(Radius)
        k = length(radius_sample)
        denominator_sample = zeros(k)
        mean_sample = zeros(k)
        y_beta_indv_sample = zeros(k, param.deg)

        #Calculate mean:
        for i in 1:k
            results = cond_density_estimation(y=radius_sample[i], y_max=param.Radius_max, y_min=param.Radius_min, x_max=param.Mass_max, x_min=param.Mass_min, deg=param.deg, w_hat=param.weights_mle, qtl=quantile, only_output_mean=true) ###what function is this??? Is it part of some CRAN package...
            mean_sample[i] = (results[1]) ###why are there parentheses here? It makes no difference...
            denominator_sample[i] = results[2]
            y_beta_indv_sample[i,:] = results[3:57]
        end
        predicted_mean = mean(mean_sample)

        #Calculate 16% and 84% quantiles:
        #Mixture of the CDF of k conditional densities
        function pbeta_conditional_density(x)
            function mix_density(j)
                deg_vec = 1:param.deg
                x_indv_cdf = map(degree -> cdf.(Beta(degree, param.deg-degree+1), (j-x_min)/(x_max-x_min)), deg_vec)
                quantile_numerator = zeros(k)
                p_beta_sample = zeros(k)
                for ii in 1:k
                    quantile_numerator[ii] = sum(param.weights_mle .* kron(x_indv_cdf, y_beta_indv_sample[ii,:]))
                    p_beta_sample[ii] = quantile_numerator[ii]/denominator_sample[ii]
                end
                return p_beta = mean(p_beta_sample)
            end

            return map(mix_density, x)
        end

        function mixture_conditional_quantile(q, x_min, x_max)
            function g(x)
                return pbeta_conditional_density(x) - q
            end
            return find_zero(g, (x_min, x_max)) ###see note above
            #return Optim.minimum(optimize(g, x_min, x_max))
        end

        predicted_quantiles = map(q -> mixture_conditional_quantile(q, param.Mass_min, param.Mass_max), qtl)
        #predicted_lower_quantile = predicted_quantiles[1]
        #predicted_upper_quantile = predicted_quantiles[2]
    end

    #Return the output:
return (predicted_mean, predicted_quantiles)::Tuple{Float64, Vector{Float64}}
end

function draw_planet_mass_from_radius_Ning2018(Radius::Float64, param::MR_param_Ning2018)
    #This function takes in a Radius (in solar radii) and draws a mass (returning in solar masses) probabilistically
    Radius_in_earths = Radius/ExoplanetsSysSim.earth_radius
    @assert param.Radius_min < log10.(Radius_in_earths) < param.Radius_max
    q::Float64 = rand()
    l_mass = predict_mass_given_radius(Radius_in_earths, param; qtl=[q])[2][1]
    return (10^l_mass)*ExoplanetsSysSim.earth_mass
end

function generate_planet_mass_from_radius_Ning2018(Radius::Float64, sim_param::SimParam)
    global MR_param
    return draw_planet_mass_from_radius_Ning2018(Radius, MR_param)
end





##### To initialize the M-R model parameters:

weights_mle = CSV.read("weights.mle.csv")[:x]
degrees = 55

N_keep = 25
indices_keep = sortperm(weights_mle, rev=true)[1:N_keep]
indices_div = 1+div.(indices_keep .- 1, degrees)
indices_rem = 1+rem.(indices_keep .- 1, degrees)
weights = weights_mle[indices_keep]
MR_param = MR_param_Ning2018(-1., 3.809597, -0.302, 1.357509, degrees, weights, indices_keep, indices_div, indices_rem) #(-1., 3.809597, -0.3, 1.357509, 55, weights_mle)





##### Examples:

#=
#Observation without measurement errors:

Radius = 5. #original scale, not log scale
predict_result = predict_mass_given_radius(Radius, MR_param)
println(predict_result)

#Observation with a measurement error:

Radius = 5. #original scale, not log scale
R_sigma = 0.1
predict_result = predict_mass_given_radius(Radius, MR_param; R_sigma=0.1)
println(predict_result)

#If want to change 16% and 84% quantiles to 5% and 95% quantiles:

Radius = 5. #original scale, not log scale
predict_result = predict_mass_given_radius(Radius, MR_param; qtl=[0.05, 0.95])
println(predict_result)
=#

#Input are posterior samples: ###currently broken because the function 'cond_density_estimation' is undefined
#=
Radius_dist = Normal(5, 0.5) #original scale, not log scale
Radius_samples = rand(Radius_dist, 100)
predict_result = predict_mass_given_radius(Radius_samples, MR_param; R_sigma=0.1, posterior_sample=true)
println(predict_result)
=#





##### For timing the function:
#=
sim_param = setup_sim_param_model()
#Radii = ones(10000)*ExoplanetsSysSim.earth_radius
Radii = (10.^(linspace(MR_param.Radius_min+0.01, MR_param.Radius_max-0.01, 10000)))*ExoplanetsSysSim.earth_radius
tic()
Masses = map(r -> generate_planet_mass_from_radius_Ning2018(r, sim_param), Radii)
t_elapsed = toc()
=#



##### For plotting and timing the model:
#=
#Pkg.add("PyPlot")
using PyPlot

N_keep_array = [3025, 19, 18]
colors = ["blue", "green", "red"]
Radii = 10.^(linspace(MR_param.Radius_min+0.01, MR_param.Radius_max-0.01, 1000))

fig, ax = subplots()
for (i,N_keep) in enumerate(N_keep_array)
    indices_keep = sortperm(weights_mle, rev=true)[1:N_keep]
    weights = weights_mle[indices_keep]
    MR_param = MR_param_Ning2018(-1., 3.809597, -0.302, 1.357509, degrees, weights, indices_keep)

    tic()
    Masses_predict = map(r -> predict_mass_given_radius(r, MR_param; qtl=[0.16, 0.84]), Radii)
    t_elapsed = toc()

    Masses_mean = 10.^[x[1] for x in Masses_predict]
    Masses_lower = 10.^[x[2][1] for x in Masses_predict]
    Masses_upper = 10.^[x[2][2] for x in Masses_predict]

    loglog(Radii, Masses_mean, color=colors[i], linewidth=2.0, linestyle="-", label="M = "*string(N_keep))
    loglog(Radii, Masses_lower, color=colors[i], linewidth=1.0, linestyle="--")
    loglog(Radii, Masses_upper, color=colors[i], linewidth=1.0, linestyle="--")
end
xlabel(L"$R_p (R_\oplus)$")
ylabel(L"$M_p (M_\oplus)$")
legend(loc="upper left")
=#





##### If we want to generate and save a pre-computed table of masses on a radius vs. quantile grid:
#=
N_keep = 3025
indices_keep = sortperm(weights_mle, rev=true)[1:N_keep]
weights = weights_mle[indices_keep]
MR_param = MR_param_Ning2018(-1., 3.809597, -0.3, 1.357509, degrees, weights, indices_keep)

N_radii = 1001
N_quantiles = 1001
log_Radii = collect(linspace(MR_param.Radius_min, MR_param.Radius_max, N_radii))
Radii = 10.^log_Radii
quantiles = collect(linspace(0., 1.0, N_quantiles))

file_name = "MRpredict_table_weights"*string(N_keep)*"_R"*string(N_radii)*"_Q"*string(N_quantiles)*".txt"
f = open(file_name, "w")
println(f, "# All masses are in log10; weights used = ", N_keep)
println(f, "# First uncommented line is the header with the column labels (first column contains the radii in log10 while the remaining columns correspond to the quantiles)")
writedlm(f, reshape(append!(["log_R"], string.(quantiles)), (1,:)), ", ")

log_Mass_table = zeros(N_radii, N_quantiles + 1)
for (i,r) in enumerate(Radii)
    log_Mass_quantiles = predict_mass_given_radius(r, MR_param; qtl=quantiles)[2]
    log_Mass_table[i,:] = append!([log_Radii[i]], log_Mass_quantiles)
    println(i, ", Radius = ", r)
end
writedlm(f, log_Mass_table, ", ")
close(f)
=#


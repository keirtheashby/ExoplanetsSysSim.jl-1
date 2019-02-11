#README
#calc_ttvfast takes a planetary system, plots TTV vs transit time, and returns an array of arrays holding the TTVs for each planet. It is slower than calc_ttvfaster but works for more systems.
#calc_ttvfaster takes a planetary system, plots TTV vs transit time, and returns an array of arrays holding the TTVs for each planet. 
#plot_ttv takes an array of TTVs, an array of transit times, and an array of periods for each planet in a system, and plots TTVs vs. transit_times.
#plot_ttv takes an array of TTVs, an array of periods, and an array of first-transit-times (t0) for each planet in a system, and creates an array of transit times before calling plot_ttvfast.
#get_ttv_sum_stats takes an array of TTVs for each planet in a system and returns the scatter, rms and iqr for all the ttvs associated with each planet.

#This code was tested outside of SysSim, to incorporate into SysSim, put the following lines in a "calc_transit_obs" function just after the for loop over each planet.
#    if num_planets(sys)>1
#      sysTTVs = calc_ttvfaster(sys) #,sysTTVs,sys_stats)		#returns an array of arrays of TTVs for each planet calculated.
#      sysTTVs = calc_ttvfast(sys) #,sysTTVs,sys_stats)               	#returns an array of arrays of TTVs for each planet calculated.
#      ttv_sum_stat = get_ttv_sum_stats(sysTTVs)			#returns the scatter, rms and iqr for all the ttvs associated with each planet.
#    end


include("TTVFaster/src/TTVFaster.jl")
include("TTVFast/src/TTVFast.jl")
#include("Stats/src/Stats.jl") required for calculating iqr but the package is currently out of date. 
using PyPlot
using StatsBase

# Integration parameters
inflag = 0
t_start=0.0
t_stop=4*365.2425
dt = 0.02
kep_obs_duration = t_stop-t_start
scatter_const = 1.4826; 

            
function calc_ttvfast(sys::Any) #sys::ExoplanetsSysSim.PlanetarySystem{ExoplanetsSysSim.Star}
#calculates transit timing variations for each planet in a system over the course of Kepler's observations. 
  nplanets = num_planets(sys)
  p = Array{Cdouble}(undef,2+nplanets*7)		#array holding information for each planet in the system
  p[1] = 4pi/365.2425^2; 				# G
  p[2] = 1.0; 						# Mstar
  num_events = 0;
  for i in 1:nplanets
    period = sys.orbit[i].P
    p[2+7*(i-1)+1] = sys.planet[i].mass 		# Planet Mass
    p[2+7*(i-1)+2] = period 	 			# Period
    p[2+7*(i-1)+3] = sys.orbit[i].ecc			# Eccentricity
    p[2+7*(i-1)+4] = sys.orbit[i].incl*180/pi		# Inclination
    p[2+7*(i-1)+5] = sys.orbit[i].asc_node*180/pi 	# Longitude of Ascending Node
    p[2+7*(i-1)+6] = sys.orbit[i].omega*180/pi		# Argument of Pericenter
    p[2+7*(i-1)+7] = sys.orbit[i].mean_anom*180/pi	# Mean Anomaly
    num_events += ceil(Integer, kep_obs_duration/p[2+7*(i-1)+2] +1); # /* large enough to fit all the transits calculated by the code*/
  end
  ttvfast_input = TTVFast.ttvfast_inputs_type(p, t_start=t_start, t_stop=t_stop, dt=dt)
  incl_rvs = false #if true, Make sure first rv_time is _after_ t_start+dt/2 or else won't get any RV outputs
  if incl_rvs
    num_rvs = 100
    rv_times = collect(range(t_start+0.501*dt, stop=t_stop, length=num_rvs))
    ttvfast_output = TTVFast.ttvfast_outputs_type(num_events ,rv_times)
  else
    ttvfast_output = TTVFast.ttvfast_outputs_type(num_events)
  end
  println(stderr, "# About to call TTVFast")
  @time TTVFast.ttvfast!(ttvfast_input,ttvfast_output)
  sysTTVs = Array{Array{Float64,1}}(undef,nplanets)   		#an array of all the TTVs for each planet
  transit_times = Array{Array{Float64}}(undef,nplanets)		#arrays for each planet holding actual transit time as measured from the beginning of kepler observation.
  ttv_per = zeros(nplanets)					#stores period values for plotting
  for k in 1:nplanets
    sysTTVs[k] = [0.0]
    transit_times[k]=[0.0]
  end
  for j in 1:num_events						#extract transit times
    this_planet = TTVFast.get_event(ttvfast_output,j).planet+1	#determine what planet we are extracting ttvs for
    this_time = TTVFast.get_event(ttvfast_output,j).time	#extract ttv time
    if this_time != 0.0						#events marked 0.0 are extras
      push!(transit_times[this_planet],this_time) 
    end
  end 
  slope = 0.0
  intercept = 0.0
  for k in 1:nplanets					
    popfirst!(transit_times[k])					#remove initializing 0.0 values 
    num_times = length(transit_times[k])
    #perform linear regression on transit_times to determine actual period (slope) and t0 (intercept)
    xmean = (num_times-1)/2
    ymean = StatsBase.mean(transit_times[k])
    numer = 0.0
    denom = 0.0
    for m in 1:num_times
      numer += ((m-1)-xmean)*(transit_times[k][m]-ymean)
      denom += ((m-1)-xmean)^2
    end
    slope = numer/denom
    ttv_per[k] = slope
    intercept = ymean-slope*xmean
    for n in 1:num_times
      this_ttv = transit_times[k][n] - slope*(n-1) - intercept	#calculate ttv time relative to expected transit time
      push!(sysTTVs[k],this_ttv)
    end
    popfirst!(sysTTVs[k]) 					#removes 0.0 value from when we initialized our array
  end
  TTVFast.free_ttvfast_outputs(ttvfast_output)			#original ttv code recommends running this 
  plot_ttv(transit_times,sysTTVs,ttv_per)			#plots relative ttvs over total observation time 
  return sysTTVs 
end


function calc_ttvfaster(sys::Any) #ExoplanetsSysSim.PlanetarySystem{ExoplanetsSysSim.Star}
#takes information from the second planet in a pair of two, and returns an array of arrays of TTVs for each planet. 
  nplanets = num_planets(sys)
  sysTTVs = Array{Array{Float64,1}}(undef,nplanets)   #an array of all the TTVs for each planet
  ttv_per = zeros(nplanets)
  ttv_t0 = zeros(nplanets)
  for p in 2:nplanets
    total_transits = 1550
    jmax = 5
    mu1 = sys.planet[p-1].mass/sys.star.mass					#use reduced mass of planet and star
    mu2 = sys.planet[p].mass/sys.star.mass
    per1 = sys.orbit[p-1].P
    per2 = sys.orbit[p].P
    t01 = rand(0:per1)								#we'll use this until we have a better function for initial transit time
    t02 = rand(0:per2)
    ecc1 = sys.orbit[p-1].ecc
    ecc2 = sys.orbit[p].ecc
    omega1 = sys.orbit[p-1].omega
    omega2 = sys.orbit[p].omega
    p1 = TTVFaster.Planet_plane_hk(mu1,per1,t01,ecc1*cos(omega1),ecc1*sin(omega1)) #organize data for each planet
    p2 = TTVFaster.Planet_plane_hk(mu2,per2,t02,ecc2*cos(omega2),ecc2*sin(omega2))
    n1 = convert(Int64,floor(total_transits/per1))				#number of transits for each planet
    n2 = convert(Int64,floor(total_transits/per2))
    time1 = collect(p1.trans0 .+ range(0,stop=n1-1,length=n1) .* p1.period)
    time2 = collect(p2.trans0 .+ range(0,stop=n2-1,length=n2) .* p2.period)
    alpha0=(p1.period/p2.period)^(2//3)
    # Initialize the computation of the Laplace coefficients:
    b0=TTVFaster.LaplaceCoefficients.initialize(jmax+1,alpha0)
    # Define arrays to hold the TTVs:
    data_array = [mu1,per1,t01,ecc1*cos(omega1),ecc1*sin(omega1),mu2,per2,t02,ecc2*cos(omega2),ecc2*sin(omega2)]
    ttv_el_type = eltype(data_array) == Float64 ? Float64 : Number
    ttv1=Array{ttv_el_type}(undef,n1)
    ttv2=Array{ttv_el_type}(undef,n2)
    # Define arrays to hold the TTV coefficients and Laplace coefficients:
    f1=Array{Float64}(undef,jmax+2,5)
    f2=Array{Float64}(undef,jmax+2,5)
    b1 = Array{Float64}(undef,jmax+2,3)               #b1 was renamed from b because b is impact parameter
    TTVFaster.compute_ttv!(jmax,p1,p2,time1,time2,ttv1,ttv2,f1,f2,b1,alpha0,b0)
    #store TTVs in sysTTVs
    if p>2
      for h in 1:n1
        sysTTVs[p-1][h] += ttv1[h] #combine ttvs of planet 1 on planet 2 and planet 3 on planet 2
      end
      sysTTVs[p] = ttv2
      ttv_per[p] = per2
      ttv_t0[p] = t02
    else
      sysTTVs[1] = ttv1
      sysTTVs[2] = ttv2
      ttv_per[1] = per1
      ttv_per[2] = per2
      ttv_t0[1] = t01
      ttv_t0[2] = t02 
    end
  end
  plot_ttvfaster(sysTTVs,ttv_per,ttv_t0)                         #plots data returned by calc_ttvfaster. Does not work with calc_ttvfast.
  return sysTTVs
end

function plot_ttv(transit_time::Array{Array{Float64,N} where N,1}, TTVs::Array{Array{Float64,1},1}, ttv_per::Array{Float64,1} )
#testing function that plots relative ttvs over total observation time.
  for i in 1:length(TTVs)
    this_per = round(Int64,ttv_per[i])
    plot(transit_time[i],TTVs[i], "-o", ms = 5, label = "period: $this_per")
    ylabel("TTV (days)", fontsize = 14)
    xlabel("time (days)", fontsize = 14)
    legend()
  end
end

function plot_ttvfaster(TTVs::Array{Array{Float64,1},1}, ttv_per::Array{Float64,1}, ttv_t0::Array{Float64,1})
  transit_time = Array{Array{Float64}}(undef,length(TTVs))
  for k in 1:length(TTVs)
    first_time = ttv_t0[1] + TTVs[k][1]	#we compute the first time separately so that we can overwrite the "undefined" initializer and then push the remaining times
    transit_time[k] = [first_time]
    for m in 2:length(TTVs[k][:])	#start at 2 since we already filled in the first values in the line above.
      this_time = ttv_t0[k] + ttv_per[k] * (m-1) + TTVs[k][m]	 #TODO: find out if there is a way to extract a observation-based period from ttvfaster code rather than using the period generated by SysSim
      push!(transit_time[k],this_time) 
    end
  end
  plot_ttv(transit_time, TTVs, ttv_per)
end


function get_ttv_sum_stats(TTVs::Array{Array{Float64,1},1})
  num_planets = length(TTVs)
  ttv_rms = zeros(num_planets)
  ttv_iqr = zeros(num_planets)
  ttv_scatter = zeros(num_planets)
  for i in 1:num_planets
    num_events = length(TTVs[i])
    ttv_mean = mean(TTVs[i])
    abs_ttvs = zeros(num_events)
    for j in 1:num_events
      abs_ttvs[j] = abs(ttv_mean - TTVs[i][j])
    end
    ttv_scatter[i] = mean(abs_ttvs)*scatter_const
    #ttv_iqr[i] = iqr(absttvs)
    ttv_rms[i] = convert(Float64,sqrt(ttv_mean^2))
  end
  ttv_sum_stat = Dict("scatter" => ttv_scatter, "rms" => ttv_rms, "iqr" => ttv_iqr)
  #println(ttv_sum_stat)
  return ttv_sum_stat
end



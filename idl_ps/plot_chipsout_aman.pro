
;; IDL code to take the output binary files from CHIPS and plot them with Bryna's plotting code. The code computes and applies the normalisations to take the units from Jy^2Hz^2sr^2 to mK^2 Mpc^3 and transposes and folds the data across the kpar = 0 line. The code will plot the output from two runs (or two pols)

; Enter details

band = 'high'    ; options: 'high', 'low'

plot_max_1D = 1.e8
plot_min_1D = 1.e2 ;1e-1

plot_max_2D = 1.e13
plot_min_2D = 1.e3

num_k_bins = 30
low_k_bin_1D = 5e-2 ;5e-3

wedge = 3.   ; [0 = no cut; 3.0 = horizon]

kperp_min_1D = 0.001
kperp_max_1D = 3.0

kpar_min_1D = 0.
kpar_max_1D = 10.

basename= '/fred/oz048/achokshi/mwa_dipole/simulations/full_sim_II/chips/AC_1125953248_full_II/'
;instring1 = 'xx_0.iter.ultralow_test2_BP2'
;instring2 = 'xx_0.iter.ultralow_test_iono_amps_SY'
; 'xx_0.iter.ultralow_test_SY_BP1'
; 'xx_0.iter.ultralow_test2_BP2'
; 'xx_0.iter.ultralow_test2_avg_BP3'
; 'xx_0.iter.ultralow_test2_fit_BP4'


plotting_mode = ''
outputstring = ''
mxval = ''
instring1 = ''
instring2 = ''
;read,basename,prompt='Enter the full path for the CHIPS output files: '
read,instring1,prompt='Enter the full extension for the first set of files: '
read,instring2,prompt='Enter the full extension for the second set of files: '
read,plottingmode,prompt='Enter plotting mode (0=screen crosspowers, 1=PNG output all files, 2=screen 1d plot, 3=no plot): '
if (plottingmode eq 1) then read,outputstring,prompt='Enter output file basename: '

if (plottingmode eq 2) then read,outputstring,prompt='Enter output file basename: '
; standard dimensions output from CHIPS

nbins = 80
Nchan = 384
Nchanall = 384
chan_width = 80.e3
deltat = 8.

if (band eq 'high') then begin

z = 6.80
DM = 5990. ; transverse comoving distance z=6.8
Npoint = 9
beamxx = read_binary('../beam_area_point_xx_srhz.dat',data_type=4)
beam_point_weight = fltarr(Npoint)
beam_point_weight = [0.,1.,1.,1.,1.,1.,1.,1.,0.]   ; assume an equal weighting of pointings
beam_point_weight = beam_point_weight/total(beam_point_weight)

beam_area = total(beam_point_weight*beamxx)
beam_area_per_chan = beam_area/Nchanall

Tsys = 220. ;K 

freq = dindgen(Nchan)*chan_width + 167.035e6

endif else begin

z = 8.34
DM = 6327. ; transverse comoving distance z=8.3
Npoint = 9
beamxx = read_binary('../beam_area_point_xx_srhz.dat',data_type=4)
beam_point_weight = fltarr(Npoint)
beam_point_weight = [0.,1.,1.,1.,1.,1.,1.,1.,0.]   ; assume an equal weighting of pointings
beam_point_weight = beam_point_weight/total(beam_point_weight)

beam_area = total(beam_point_weight*beamxx)*(182./152.)^2
beam_area_per_chan = beam_area/Nchanall
freq = dindgen(Nchan)*chan_width + 138.875e6

Tsys = 220.*(182./152.)^2.6 ;K 

endelse


umax = 300.
Neta = Nchan/2+1
Neta = Nchan+1

eta_actual = dindgen(Nchan)-Nchan/2.+1.
loc = where(abs(eta_actual) lt 0.01)

tvlct,[0,255,0,0],[0,0,255,0],[0,0,0,255]
device,decomposed=0

lperp = dindgen(nbins)*umax*1.1/float(nbins)*2.*!pi
lperp[0] = lperp[1]/2.

; **** read-in data ****

crosspower = read_binary(basename+'crosspower_'+instring1+'.dat',data_type=4,data_dims=[Nchan,nbins])
totpower = read_binary(basename+'totpower_'+instring1+'.dat',data_type=4,data_dims=[Nchan,nbins])
flagpower = read_binary(basename+'flagpower_'+instring1+'.dat',data_type=4,data_dims=[Nchan,nbins])
residpower = read_binary(basename+'residpower_'+instring1+'.dat',data_type=4,data_dims=[Nchan,nbins])
weights = read_binary(basename+'outputweights_'+instring1+'.dat',data_type=4,data_dims=[Nchan,nbins])
;fg_power = read_binary(basename+'fg_power_'+instring1+'.dat',data_type=4,data_dims=[Nchan,nbins])
fg_num = read_binary(basename+'fg_num_'+instring1+'.dat',data_type=3,data_dims=[Nchan,nbins])


crosspower_2 = read_binary(basename+'crosspower_'+instring2+'.dat',data_type=4,data_dims=[Nchan,nbins])
totpower_2 = read_binary(basename+'totpower_'+instring2+'.dat',data_type=4,data_dims=[Nchan,nbins])
flagpower_2 = read_binary(basename+'flagpower_'+instring2+'.dat',data_type=4,data_dims=[Nchan,nbins])
residpower_2 = read_binary(basename+'residpower_'+instring2+'.dat',data_type=4,data_dims=[Nchan,nbins])
weights_2 = read_binary(basename+'outputweights_'+instring2+'.dat',data_type=4,data_dims=[Nchan,nbins])
;fg_power_2 = read_binary(basename+'fg_power_'+instring2+'.dat',data_type=4,data_dims=[Nchan,nbins])
fg_num_2 = read_binary(basename+'fg_num_'+instring2+'.dat',data_type=3,data_dims=[Nchan,nbins])

fg_power = crosspower*0.
fg_power_2 = crosspower*0.


;instring1 = 'EW good_200'
;instring2 = 'EW bad_200'

; **** start transposing and folding data ****


crosspower = transpose(crosspower)
totpower = transpose(totpower)
flagpower = transpose(flagpower)
residpower = transpose(residpower)
weights = transpose(weights)
fg_power = transpose(fg_power)
fg_num = transpose(fg_num)

crosspower_2 = transpose(crosspower_2)
totpower_2 = transpose(totpower_2)
flagpower_2 = transpose(flagpower_2)
residpower_2 = transpose(residpower_2)
weights_2 = transpose(weights_2)
fg_power_2 = transpose(fg_power_2)
fg_num_2 = transpose(fg_num_2)

crosspower_temp = fltarr(nbins,Nchan/2)
totpower_temp = fltarr(nbins,Nchan/2)
flagpower_temp = fltarr(nbins,Nchan/2)
residpower_temp = fltarr(nbins,Nchan/2)
;residpowerimag_temp = fltarr(nbins,Nchan/2)
weights_temp = fltarr(nbins,Nchan/2)
fg_power_temp = fltarr(nbins,Nchan/2)
fg_num_temp = fltarr(nbins,Nchan/2)

fold=0
if (fold eq 1) then begin

; fold data

for i=0,nbins-1 do begin

	crosspower_temp[i,0] = crosspower[i,loc]
	flagpower_temp[i,0] = flagpower[i,loc]
	residpower_temp[i,0] = residpower[i,loc]
	totpower_temp[i,0] = totpower[i,loc]
	weights_temp[i,0] = weights[i,loc]
    	fg_power_temp[i,0] = fg_power[i,loc]
        fg_num_temp[i,0] = fg_num[i,loc]


    for j=1,Nchan/2-1 do begin

    crosspower_temp[i,j] = 0.5*(crosspower[i,j+loc] + crosspower[i,loc-j])
    totpower_temp[i,j] = 0.5*(totpower[i,j+loc] + totpower[i,loc-j])
    flagpower_temp[i,j] = 0.5*(flagpower[i,j+loc] + flagpower[i,loc-j])
    residpower_temp[i,j] = 0.5*(residpower[i,j+loc] + residpower[i,loc-j])
    weights_temp[i,j] = 0.5*(weights[i,j+loc] + weights[i,loc-j])
    fg_power_temp[i,j] = 0.5*(fg_power[i,j+loc] + fg_power[i,loc-j])
    fg_num_temp[i,j] = 0.5*(fg_num[i,j+loc] + fg_num[i,loc-j])


    endfor
endfor



crosspower = crosspower_temp
residpower = residpower_temp
totpower = totpower_temp
weights = weights_temp
flagpower = flagpower_temp
fg_power = fg_power_temp
fg_num = fg_num_temp


; String 2

; fold data

for i=0,nbins-1 do begin

	crosspower_temp[i,0] = crosspower_2[i,loc]
	flagpower_temp[i,0] = flagpower_2[i,loc]
	residpower_temp[i,0] = residpower_2[i,loc]
;	residpowerimag_temp[i,0] = residpowerimag_2[i,loc]
	totpower_temp[i,0] = totpower_2[i,loc]
	weights_temp[i,0] = weights_2[i,loc]
    	fg_power_temp[i,0] = fg_power_2[i,loc]
        fg_num_temp[i,0] = fg_num_2[i,loc]

    for j=1,Nchan/2-1 do begin

    crosspower_temp[i,j] = 0.5*(crosspower_2[i,j+loc] + crosspower_2[i,loc-j])
    totpower_temp[i,j] = 0.5*(totpower_2[i,j+loc] + totpower_2[i,loc-j])
    flagpower_temp[i,j] = 0.5*(flagpower_2[i,j+loc] + flagpower_2[i,loc-j])
    residpower_temp[i,j] = 0.5*(residpower_2[i,j+loc] + residpower_2[i,loc-j])
;    residpowerimag_temp[i,j] = 0.5*(residpowerimag_2[i,j+loc] + residpowerimag_2[i,loc-j])
    weights_temp[i,j] = 0.5*(weights_2[i,j+loc] + weights_2[i,loc-j])
    fg_power_temp[i,j] = 0.5*(fg_power_2[i,j+loc] + fg_power_2[i,loc-j])
    fg_num_temp[i,j] = 0.5*(fg_num_2[i,j+loc] + fg_num_2[i,loc-j])

    endfor
endfor



crosspower_2 = crosspower_temp
residpower_2 = residpower_temp
totpower_2 = totpower_temp
weights_2 = weights_temp
flagpower_2 = flagpower_temp
fg_power_2 = fg_power_temp
fg_num_2 = fg_num_temp


endif else begin

for i=0,nbins-1 do begin

	crosspower_temp[i,0] = crosspower[i,loc]
	flagpower_temp[i,0] = flagpower[i,loc]
	residpower_temp[i,0] = residpower[i,loc]
	totpower_temp[i,0] = totpower[i,loc]
	weights_temp[i,0] = weights[i,loc]
    	fg_power_temp[i,0] = fg_power[i,loc]
        fg_num_temp[i,0] = fg_num[i,loc]


    for j=1,Nchan/2-1 do begin

    crosspower_temp[i,j] = (crosspower[i,j+loc])
    totpower_temp[i,j] = (totpower[i,j+loc])
    flagpower_temp[i,j] = (flagpower[i,j+loc])
    residpower_temp[i,j] = (residpower[i,j+loc])
    weights_temp[i,j] = (weights[i,j+loc])
    fg_power_temp[i,j] = (fg_power[i,j+loc])
    fg_num_temp[i,j] = (fg_num[i,j+loc])


    endfor
endfor

crosspower = crosspower_temp
residpower = residpower_temp
totpower = totpower_temp
weights = weights_temp
flagpower = flagpower_temp
fg_power = fg_power_temp
fg_num = fg_num_temp


; FHD

; fold data

for i=0,nbins-1 do begin

	crosspower_temp[i,0] = crosspower_2[i,loc]
	flagpower_temp[i,0] = flagpower_2[i,loc]
	residpower_temp[i,0] = residpower_2[i,loc]
;	residpowerimag_temp[i,0] = residpowerimag_2[i,loc]
	totpower_temp[i,0] = totpower_2[i,loc]
	weights_temp[i,0] = weights_2[i,loc]
    	fg_power_temp[i,0] = fg_power_2[i,loc]
        fg_num_temp[i,0] = fg_num_2[i,loc]

    for j=1,Nchan/2-1 do begin

    crosspower_temp[i,j] = (crosspower_2[i,j+loc])
    totpower_temp[i,j] = (totpower_2[i,j+loc])
    flagpower_temp[i,j] = (flagpower_2[i,j+loc])
    residpower_temp[i,j] = (residpower_2[i,j+loc])
;    residpowerimag_temp[i,j] = (residpowerimag_2[i,j+loc])
    weights_temp[i,j] = (weights_2[i,j+loc])
    fg_power_temp[i,j] = (fg_power_2[i,j+loc])
    fg_num_temp[i,j] = (fg_num_2[i,j+loc])

    endfor
endfor



crosspower_2 = crosspower_temp
residpower_2 = residpower_temp
;residpowerimag_2 = residpowerimag_temp
totpower_2 = totpower_temp
weights_2 = weights_temp
flagpower_2 = flagpower_temp
fg_power_2 = fg_power_temp
fg_num_2 = fg_num_temp

endelse

; **** end transposing and folding data ****


; *****************************

; Define parameters

; define cosmological units and conversions

;deltat = 8.
hubble = 0.704*100.  ; km/s/Mpc
c = 2.995e5   ; km/s
omega_matter = 0.272
omega_baryon = 0.046
omega_lambda = 0.7
n = 0.5  ; scale-invariant spectrum
deltah = 1.94e-5*omega_matter^(-0.785-0.05*alog(omega_matter))*0.93
bw = Nchan*chan_width   ; Hz
f21 = c*1.e3/0.21
Ez = sqrt(omega_matter*(1.+z)^3 + omega_lambda)

boltz = 1.38e-23
D = 4.6  ;m
lambda = 2.997e8/freq   ;m
; beam area calculated from int (beam^2) for each coarse channel

beam_area_steradian = beam_area_per_chan/chan_width

obs_volume = mean(beam_area_steradian)*chan_width*Nchan    ; sr. Hz

jy2_Nchan__to__jy2hz2 = chan_width^2*Nchanall
jy2__to__K2sr2 = (mean(lambda^2)/(2.*boltz*1.e26))^2

mpc_conversion = DM^2*2.995e5*(1.+z)^2/100./f21/Ez    ; Mpc^3/sr.Hz

normalisation = jy2_Nchan__to__jy2hz2*jy2__to__K2sr2/obs_volume*1.e6*mpc_conversion

M = Nchan/2

;normm = 177.^2/64.;*sqrt(2.)
normm=1.

normalisation_w = normalisation

expected_noise = 2.*boltz*Tsys*1.e26/D^2/sqrt(bw/Nchan*deltat)   ; Jy
expected_noise = (expected_noise)^2   ; square ->Jy^2

factor = c*(1.+z)^2/(2.*!pi*100.*f21*Ez)

; parametrize in h units - (/Mpc -> h /Mpc)
hfactor = hubble/100. ; Hubble constant used in calcs.

; ************************************************

pratio = crosspower*0.0
flagratio = crosspower*0.0
ptotratio = crosspower*0.0
residratio = crosspower*0.0
kkk = where(weights ne 0.)
pratio[kkk] = crosspower[kkk]/weights[kkk]*normalisation
flagratio[kkk] = flagpower[kkk]/weights[kkk]*normalisation*expected_noise
ptotratio[kkk] = totpower[kkk]/weights[kkk]*normalisation
residratio[kkk] = residpower[kkk]/weights[kkk]*normalisation

weight_data = weights/(normalisation_w)^2*(440./Tsys)^4;*(152./182.)^5.2;/16.;/4./4.

;normm=1.
crosspower = pratio/normm^2/4.
residpower = residratio/normm^2/4.
totpower = ptotratio/normm^2/4.
weights = weight_data*normm^4*16.
flagpower = flagratio/normm^2/4.

; ********

pratio = crosspower*0.0
flagratio = crosspower*0.0
ptotratio = crosspower*0.0
residratio = crosspower*0.0
;residimagratio = crosspower*0.0
kkk = where(weights_2 ne 0.)
pratio[kkk] = crosspower_2[kkk]/weights_2[kkk]*normalisation
flagratio[kkk] = flagpower_2[kkk]/weights_2[kkk]*normalisation*expected_noise
ptotratio[kkk] = totpower_2[kkk]/weights_2[kkk]*normalisation
residratio[kkk] = residpower_2[kkk]/weights_2[kkk]*normalisation

weight_data = weights_2/(normalisation_w)^2*(440./Tsys)^4;*(152./182.)^5.2;/16.;/4./4.

;normm=1.
crosspower_2 = pratio/normm^2/4.
residpower_2 = residratio/normm^2/4.
;residpowerimag_2 = residimagratio
totpower_2 = ptotratio/normm^2/4.
weights_2 = weight_data*normm^4*16.
flagpower_2 = flagratio/normm^2/4.

weight_data = 0.
pratio=0.

; setup ratios and differences between the two inputs

kgood = where(crosspower_2 ne 0.)

ratio_cross = crosspower*0.

ratio_cross[kgood] =  crosspower[kgood]/crosspower_2[kgood]
diff_cross =  crosspower - crosspower_2

crosspower = (crosspower)

; ********************************************************************

eta = findgen(Neta)-0.5   ; Using Bryna's definition for zeroth bin

kpa = eta/bw/factor
kper = lperp/DM

conv_factor = 1./(1./DM*2.*!pi)
ns_conv = eta/bw*1.e9;/hfactor


; 1.06 here gives the horizon line for the maximum pointing
wedgeedge = 1.06*DM/2./!pi/factor/freq[0];/hfactor


mx = plot_max_2D
mn = plot_min_2D

Nchancut=Nchan

delay_params = [ns_conv[2]-ns_conv[1] , ns_conv[Nchancut/2-1]]
kpar_bin = kpa[2]-kpa[1]

if (plottingmode eq 3) then stop

if (plottingmode eq 0) then begin

struc = {ncol:2,nrow:1,ordering:'col'}

kpower_2d_plots,kperp_edges=kper[2:nbins-10],kpar_edges=kpa[0:Nchancut/2.],/hinv,hubble_param=0.7,power=(crosspower_2[2:nbins-11,0:Nchancut/2.-1]),/plot_wedge_line,wedge_amp=[1.,wedgeedge],title_prefix=instring2,/baseline_axis,kperp_lambda_conv=conv_factor,/delay_axis,delay_params=delay_params,data_range=[mn,mx],start_multi_params=struc,weights=weights_2[2:nbins-11,0:Nchancut/2-1],kpar_bin=kpar_bin;,/png,plotfile='XX_power.png';,/plot_noise,noise_meas=(residpower_2[2:nbins-11,0:Nchancut/2]);,/snr,noise_expval=abs(flagpower[2:nbins-11,0:Nchancut/2-1])]

kpower_2d_plots,kperp_edges=kper[2:nbins-10],kpar_edges=kpa[0:Nchancut/2.],/hinv,hubble_param=0.7,power=(crosspower[2:nbins-11,0:Nchancut/2.-1]),/plot_wedge_line,wedge_amp=[1.,wedgeedge],title_prefix=instring1,/baseline_axis,kperp_lambda_conv=conv_factor,/delay_axis,delay_params=delay_params,data_range=[mn,mx],multi_pos=[0.5,0.,1.0,1.0],weights=weights[2:nbins-11,0:Nchancut/2-1],kpar_bin=kpar_bin;,plotfile='XX_power.png',/png;,noise_meas=abs(residpower[1:nbins-11,0:Nchancut/2-1]);,noise_expval=abs(flagpower_2[1:nbins-11,0:Nchancut/2-1]),kpar_bin=kpar_bin

stop
endif

; ****************************
; ****************************

if (plottingmode eq 4) then begin


struc = {ncol:2,nrow:1,ordering:'col'}

mxdiff = max(abs(diff_cross))*10.

kpower_2d_plots,kperp_edges=kper[2:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=0.7,power=ratio_cross[2:nbins-10,0:Nchancut/2-1],noise_expval=abs(flagpower[2:nbins-10,0:Nchancut/2-1]),/plot_wedge_line,wedge_amp=[1.,wedgeedge],window_num=0,full_title=instring1+' Ratio',data_range=[0.1,10.],start_multi_params=struc,kpar_bin=kpar_bin,/delay_axis,delay_params=delay_params,/no_units
kpower_2d_plots,kperp_edges=kper[2:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=0.7,power=diff_cross[2:nbins-10,0:Nchancut/2-1],/plot_wedge_line,wedge_amp=[1.,wedgeedge],window_num=0,full_title=instring1+' Diff',data_range=[-mxdiff,mxdiff],multi_pos=[0.5,0.,1.0,1.0],kpar_bin=kpar_bin,/delay_axis,delay_params=delay_params,weights=weights[2:nbins-10,0:Nchancut/2-1],color_profile='sym_log'

endif



if (plottingmode eq 1) then begin

struc = {ncol:2,nrow:1,ordering:'col'}
output = 'plots_'+outputstring+'.png'

kpower_2d_plots,kperp_edges=kper[1:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=0.7,power=(crosspower[1:nbins-10,0:Nchancut/2-1]),/plot_wedge_line,wedge_amp=[1.,wedgeedge],window_num=0,title_prefix=instring1,start_multi_params=struc,plotfile=output,/png,data_range=[plot_min_2D,plot_max_2D],/baseline_axis,kperp_lambda_conv=conv_factor,/delay_axis,delay_params=delay_params,kpar_bin=kpar_bin

kpower_2d_plots,kperp_edges=kper[1:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=0.7,power=(crosspower_2[1:nbins-10,0:Nchancut/2-1]),/plot_wedge_line,wedge_amp=[1.,wedgeedge],window_num=0,multi_pos=[0.5,0.,1.0,1.0],title_prefix=instring2,plotfile=output,/png,data_range=[plot_min_2D,plot_max_2D],/baseline_axis,kperp_lambda_conv=conv_factor,/delay_axis,delay_params=delay_params,kpar_bin=kpar_bin

; ********************************

;outputstring = 'EW_8sec_48sec'
;kgood = where(crosspoweryy_96 gt 0.)
;ratio_cross[kgood] =  crosspoweryy_192[kgood]/crosspoweryy_96[kgood]
;diff_cross =  crosspoweryy_192-crosspoweryy_96


struc = {ncol:2,nrow:1,ordering:'col'}
output = 'plots_'+outputstring+'_ratiodiff.png'

mxdiff = max(abs(diff_cross))*10.

kpower_2d_plots,kperp_edges=kper[2:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=0.7,power=ratio_cross[2:nbins-10,0:Nchancut/2-1],noise_expval=abs(flagpower[2:nbins-10,0:Nchancut/2-1]),/plot_wedge_line,wedge_amp=[1.,wedgeedge],window_num=0,full_title=instring1+' Ratio',plotfile=output,/png,data_range=[0.1,10.],start_multi_params=struc,kpar_bin=kpar_bin,/delay_axis,delay_params=delay_params,/no_units
kpower_2d_plots,kperp_edges=kper[2:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=0.7,power=diff_cross[2:nbins-10,0:Nchancut/2-1],/plot_wedge_line,wedge_amp=[1.,wedgeedge],window_num=0,full_title=instring1+' Diff',plotfile=output,/png,data_range=[-mxdiff,mxdiff],multi_pos=[0.5,0.,1.0,1.0],kpar_bin=kpar_bin,/delay_axis,delay_params=delay_params,weights=weights[2:nbins-10,0:Nchancut/2-1],color_profile='sym_log'

; ********************************


struc = {ncol:3,nrow:2,ordering:'col'}
output = 'plots_'+outputstring+'_noise.png'

kpower_2d_plots,kperp_edges=kper[2:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=0.7,power=crosspower[2:nbins-10,0:Nchancut/2-1],/plot_wedge_line,wedge_amp=[1.,wedgeedge],window_num=0,title_prefix=instring1,/delay_axis,delay_params=delay_params,start_multi_params=struc,weights=weights[2:nbins-10,0:Nchancut/2-1],/plot_sigma,plotfile=output,/png,data_range=[1.e4,1.e13],kpar_bin=kpar_bin
kpower_2d_plots,kperp_edges=kper[2:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=0.7,power=(crosspower[2:nbins-10,0:Nchancut/2-1]),/plot_wedge_line,wedge_amp=[1.,wedgeedge],window_num=0,noise_expval=abs(flagpower[2:nbins-10,0:Nchancut/2-1]),/plot_exp_noise,title_prefix=instring1,multi_pos=[0.33,0.5,0.67,1],plotfile=output,/png,data_range=[1.e4,1.e13],kpar_bin=kpar_bin,/delay_axis,delay_params=delay_params
kpower_2d_plots,kperp_edges=kper[2:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=0.7,power=crosspower[2:nbins-10,0:Nchancut/2-1],/plot_wedge_line,wedge_amp=[1.,wedgeedge],window_num=0,title_prefix=instring1,multi_pos=[0.67,0.5,1.0,1.0],noise_meas=(residpower[2:nbins-10,0:Nchancut/2-1]),/plot_noise,plotfile=output,/png,data_range=[1.e4,1.e13],kpar_bin=kpar_bin,/delay_axis,delay_params=delay_params

kpower_2d_plots,kperp_edges=kper[2:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=0.7,power=crosspower_2[2:nbins-10,0:Nchancut/2-1],/plot_wedge_line,wedge_amp=[1.,wedgeedge],window_num=0,title_prefix=instring2,multi_pos=[0.,0.,0.33,0.5],weights=weights_2[2:nbins-10,0:Nchancut/2-1],/plot_sigma,plotfile=output,/png,data_range=[1.e4,1.e13],kpar_bin=kpar_bin,/delay_axis,delay_params=delay_params
kpower_2d_plots,kperp_edges=kper[2:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=0.7,power=(crosspower_2[2:nbins-10,0:Nchancut/2-1]),/plot_wedge_line,wedge_amp=[1.,wedgeedge],window_num=0,noise_expval=abs(flagpower_2[2:nbins-10,0:Nchancut/2-1]),/plot_exp_noise,title_prefix=instring2,multi_pos=[0.33,0.,0.67,0.5],plotfile=output,/png,data_range=[1.e4,1.e13],kpar_bin=kpar_bin,/delay_axis,delay_params=delay_params
kpower_2d_plots,kperp_edges=kper[2:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=0.7,power=crosspower_2[2:nbins-10,0:Nchancut/2-1],/plot_wedge_line,wedge_amp=[1.,wedgeedge],window_num=0,title_prefix=instring2,multi_pos=[0.67,0.,1.0,0.5],noise_meas=(residpower_2[2:nbins-10,0:Nchancut/2-1]),/plot_noise,plotfile=output,/png,data_range=[1.e4,1.e13],kpar_bin=kpar_bin,/delay_axis,delay_params=delay_params

; ****************************
struc = {ncol:2,nrow:1,ordering:'col'}
output = 'plots_'+outputstring+'_errors.png'

kpower_2d_plots,kperp_edges=kper[2:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=0.7,power=(crosspower[2:nbins-10,0:Nchancut/2-1]),/plot_wedge_line,wedge_amp=[1.,wedgeedge],window_num=0,title_prefix=instring1,start_multi_params=struc,plotfile=output,/png,data_range=[1.e4,1.e10],/baseline_axis,kperp_lambda_conv=conv_factor,/delay_axis,delay_params=delay_params,/plot_sigma,weights=weights[2:nbins-10,0:Nchancut/2-1],kpar_bin=kpar_bin

kpower_2d_plots,kperp_edges=kper[2:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=0.7,power=(crosspower_2[2:nbins-10,0:Nchancut/2-1]),/plot_wedge_line,wedge_amp=[1.,wedgeedge],window_num=0,multi_pos=[0.5,0.,1.0,1.0],title_prefix=instring2,plotfile=output,/png,data_range=[1.e4,1.e10],/baseline_axis,kperp_lambda_conv=conv_factor,/delay_axis,delay_params=delay_params,/plot_sigma,weights=weights_2[2:nbins-10,0:Nchancut/2-1],kpar_bin=kpar_bin


; ********************

struc = {ncol:2,nrow:2,ordering:'col'}
output = 'plots_'+outputstring+'_snr.png'

kpower_2d_plots,kperp_edges=kper[2:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=0.7,power=crosspower[2:nbins-10,0:Nchancut/2-1],noise_expval=abs(flagpower[2:nbins-10,0:Nchancut/2-1]),noise_meas=abs(residpower[2:nbins-10,0:Nchancut/2-1]),/nnr,/plot_wedge_line,wedge_amp=[1.,wedgeedge],window_num=0,title_prefix=instring1,plotfile=output,/png,data_range=[0.01,100.],start_multi_params=struc,kpar_bin=kpar_bin,/delay_axis,delay_params=delay_params
kpower_2d_plots,kperp_edges=kper[2:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=0.7,power=crosspower[2:nbins-10,0:Nchancut/2-1],noise_expval=abs(flagpower[2:nbins-10,0:Nchancut/2-1]),/snr,/plot_wedge_line,wedge_amp=[1.,wedgeedge],window_num=0,title_prefix=instring1,plotfile=output,/png,data_range=[1.e-2,1.e6],multi_pos=[0.,0.,0.5,0.5],kpar_bin=kpar_bin,/delay_axis,delay_params=delay_params,weights=weights[2:nbins-10,0:Nchancut/2-1]


kpower_2d_plots,kperp_edges=kper[2:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=0.7,power=crosspower_2[2:nbins-10,0:Nchancut/2-1],noise_expval=abs(flagpower_2[2:nbins-10,0:Nchancut/2-1]),/snr,/plot_wedge_line,wedge_amp=[1.,wedgeedge],window_num=0,title_prefix=instring2,plotfile=output,/png,multi_pos=[0.5,0.,1,0.5],data_range=[1.e-2,1.e6],kpar_bin=kpar_bin,/delay_axis,delay_params=delay_params,weights=weights_2[2:nbins-10,0:Nchancut/2-1]
kpower_2d_plots,kperp_edges=kper[2:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=0.7,power=crosspower_2[2:nbins-10,0:Nchancut/2-1],noise_expval=abs(flagpower_2[2:nbins-10,0:Nchancut/2-1]),noise_meas=abs(residpower_2[2:nbins-10,0:Nchancut/2-1]),/nnr,/plot_wedge_line,wedge_amp=[1.,wedgeedge],window_num=0,title_prefix=instring2,multi_pos=[0.5,0.5,1,1],plotfile=output,/png,data_range=[0.01,100.],kpar_bin=kpar_bin,/delay_axis,delay_params=delay_params


endif

; ******************************************************************************************

if (plottingmode eq 2) then begin

; form 1D power uncertainty

P_k_spher = read_binary("P_k_mK2_8mhz.dat",data_type=5,data_dims=[41])
uu = read_binary("l.dat",data_type=5)/2./!pi   ; u value

Neta = Nchan/2.
Netaa = 50.

kmax = sqrt(kper[nbins-1]^2 + kpa[Neta-1]^2)

;kmax = 3.49

ptot = dblarr(Netaa)
psig = ptot*0.
pmeas = ptot*0.
ptotfg=ptot*0.
pmeasfg=pmeas*0.
num = intarr(Netaa)
numedge = fltarr(Netaa)
valedge = fltarr(Netaa)
;print,num
mask = weights*0.0

ktot_bins = (dindgen(Netaa+1)/float(Netaa))^3.*kmax+0.05
ktot_bins = (dindgen(Netaa+1)/float(Netaa))^1.55*kmax+0.03

ktot_bins = (dindgen(Netaa+1)/float(Netaa))^1.4*kmax+low_k_bin_1D

noise_obs_fg = sqrt(weights_2);/16.
noise_obs = sqrt(weights);/16.
;noise_obs = 1./residpower/16.

power_sims = ktot_bins^(-1)

wedge = 0.

edge = dindgen(nbins)/nbins*5.    ; 2 sigma correlation length is 5 kpar bins

for i = 1 , nbins-1 do begin
   for j = 0 , Neta-1 do begin

;if (noise_obs[i,j] ne 0.) then begin


      for k = 0 , Netaa-1 do begin
         if ((sqrt(kper[i]^2 + kpa[j]^2) ge ktot_bins[k])and(sqrt(kper[i]^2 + kpa[j]^2) lt ktot_bins[k+1])) then begin
if ((kpa[j] gt kper[i]*wedge)and(kper[i] lt kperp_max_1D)and(kper[i] ge kperp_min_1D)and(kpa[j] gt kpar_min_1D)) then begin

        ptot[k] = ptot[k] + (noise_obs[i,j])^2;*maskk[i,j]
		ptotfg[k] = ptotfg[k] + (noise_obs_fg[i,j])^2
        ;pmeas[k] = pmeas[k] + (totpower[i,j])*(noise_obs[i,j])^2;*maskk[i,j]
	;pmeasfg[k] = pmeasfg[k] + (totpower_2[i,j])*(noise_obs_fg[i,j])^2
	pmeas[k] = pmeas[k] + (crosspower[i,j])*(noise_obs[i,j])^2;*maskk[i,j]
	pmeasfg[k] = pmeasfg[k] + (crosspower_2[i,j])*(noise_obs_fg[i,j])^2
    num[k] = num[k] + float(fg_num[i,j])/Nchan
                   numedge[k] = numedge[k]+1.
                   valedge[k] = valedge[k]+edge[i]
		mask[i,j] = k
       endif
         endif
      endfor
;endif

   endfor
endfor

for k = 0,Netaa-1 do psig[k] = interpol(P_k_spher,uu*2.*!pi/9144./hfactor,ktot_bins[k])

; find the average value of the +ve x error bar for the cells that enter each ktot_bin

xerrhi = valedge*0.
fin=where(numedge gt 0.)
xerrhi[fin] = valedge[fin]/numedge[fin]*(kpa[5]-kpa[4])*2.

xerrlo = dblarr(Netaa)+(kpa[5]-kpa[4])    ; for 2 sigma contained within that cell

; ptot here is the uncertainty on the dimensionless power in mK^2

pmeas = pmeas/ptot
pmeasfg=pmeasfg/ptotfg


pmeas = pmeas*ktot_bins^3/2./!pi^2
pmeasfg = pmeasfg*ktot_bins^3/2./!pi^2

kernel = 1.

ptot = (ktot_bins^3/sqrt(ptot)/2./!pi^2)*kernel
ptotfg = (ktot_bins^3/sqrt(ptotfg)/2./!pi^2)*kernel

output = 'plots_'+outputstring+'_1D.png'

window,1,xsize=700,ysize=700,retain=2

;write_png,output,tvrd(/true)

;window,1,xsize=700,ysize=700,retain=2

;thisDevice=!D.NAME
;set_plot,'ps'
;device,filename=output,xsize=30,ysize=16,/encapsulated,bits=8
;!p.font=0
;device,set_font="Times-Roman",font_size=12
;!p.thick=4 & !x.thick=4 & !y.thick=4

plot,[0,0],[0,0],charsize=1.8,charthick=1.75,thick=1.75,color=0,background=-1,title='1D Power - !7D!3!E2!N(k) = (k!E3!NP(k)/2!7p!3!E2!N) (mK!E2!N)', $
;xrange=[0.07,4],ystyle=1,yrange=[1.e5,1.e11],xtitle='log!D10!N k (h Mpc!E-1!N)',ytitle='log!D10!N !7D!3!E2!N(k) (mK!E2!N)',/ylog,/xlog,xstyle=1
;xrange=[low_k_bin_1D,4],ystyle=1,yrange=[plot_min_1D,plot_max_1D],xtitle='log!D10!N k (h Mpc!E-1!N)',ytitle='log!D10!N !7D!3!E2!N(k) (mK!E2!N)',/ylog,/xlog,xstyle=1
xrange=[low_k_bin_1D,4],ystyle=1,yrange=[plot_min_1D,plot_max_1D],xtitle='log!D10!N k (h Mpc!E-1!N)',ytitle='log!D10!N !7D!3!E2!N(k) (mK!E2!N)',/ylog,/xlog,xstyle=1
oplot,ktot_bins[1:Netaa-1],pmeas[1:Netaa-1],color=2,psym=1,thick=2,linestyle=6

oplot,ktot_bins[1:Netaa-1],pmeasfg[1:Netaa-1],color=2,thick=2,psym=1,linestyle=6

;WRITE_CSV, '1D_text_values.csv', ktot_bins[1:Netaa-1], pmeas[1:Netaa-1],pmeasfg[1:Netaa-1]

oploterror,ktot_bins[1:Netaa-1],pmeas[1:Netaa-1],xerrhi[1:Netaa-1],ptot[1:Netaa-1],errcolor=3,/hibar
oploterror,ktot_bins[1:Netaa-1],pmeas[1:Netaa-1],xerrlo[1:Netaa-1],ptot[1:Netaa-1],errcolor=3,/lobar
oplot,ktot_bins[1:Netaa-1],pmeas[1:Netaa-1],color=3,psym=1,thick=2,linestyle=6

oploterror,ktot_bins[1:Netaa-1],pmeasfg[1:Netaa-1],xerrhi[1:Netaa-1],ptotfg[1:Netaa-1],errcolor=1,/hibar
oploterror,ktot_bins[1:Netaa-1],pmeasfg[1:Netaa-1],xerrlo[1:Netaa-1],ptotfg[1:Netaa-1],errcolor=1,/lobar
oplot,ktot_bins[1:Netaa-1],pmeasfg[1:Netaa-1],color=1,thick=2,psym=1,linestyle=6
;legend,/left_legend,charsize=1.5,[instring1,instring2],color=[3,1],textcolors=[3,1],thick=2.25,psym=[1,1]

oplot,ktot_bins,ptot,color=2,thick=2
oplot,ktot_bins,ptotfg,color=2,thick=2

;device,/close
;set_x


write_png,output,tvrd(/true)

stop

struc = {ncol:1,nrow:1,ordering:'col'}
;output = 'plots_'+outputstring+'_1D.png'
output = 'plots_1D.png'

kpower_2d_plots,kperp_edges=kper[2:nbins-3],kpar_edges=kpa[0:Nchancut/2],/hinv,hubble_param=0.7,power=(mask[2:nbins-10,0:Nchancut/2-1]),/plot_wedge_line,wedge_amp=[1.,wedgeedge],window_num=0,title_prefix=instring1,start_multi_params=struc,plotfile=output,/png,data_range=[mn,mx],/baseline_axis,kperp_lambda_conv=conv_factor,/delay_axis,delay_params=delay_params,kpar_bin=kpar_bin


endif

end

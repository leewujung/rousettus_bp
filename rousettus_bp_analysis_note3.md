# Rousettus beampattern processing and modeling ---  2017


## 2017/01/16
### Move BEM results off from server
Move most BEM calculation results off from the gal-e server. Only keep the ones used in the paper:

- **20160820_bullethead-sc-colony-0.5mm**: bullethead scaled to match dimension of Ra-colony, no rotation around y-axis
- **20160829_Ra-colony-noear-0.5mm**: Ra-colony bat head *without* ears (rotated to have mouth more parallel to x-axis)
- <span style="background-color:yellow">**20160917_Ra-colony-rotear-0.5mm**</span>: Ra-colony bat head *with ears rotated slightly downward* (rotated to have mouth more parallel to x-axis) --> this slight ear rotation eliminated the 25 kHz anomaly --> this is the final mesh used in Fig. 4 in the paper

 

## 2017/03/06
### Reprocessed Rousettus bp files
* Ben found out that the mic beampattern compensation was not done correctly: the mic_bp file has the beampattern values in degrees but in `bp_proc` the angles are in radians. Wu-Jung later found out that in `bp_proc` the function `bat2mic` was problematic that the vector `mic_vec_norm` was not normalized before used to calculate the polar angle. These are corrected as of 2017/03/01 and changed were pushed up to GitHub repo.
* Wu-Jung reprocessed Rousettus bp files on 2017/03/06. Results are stored in the following folders:
	* `proc_output_rousettus_new_raw`: raw results after `bp_proc`
	* `proc_output_rousettus_new_isgood call`: results after applying `chk_isgoodcall`
	* `proc_output_rousettus_new_checked`: results after copying the good calls from what Ben had checked before (`proc_output_rousettus_checked`). The copying was done using code `copy_isgoodcall`



## 2017/09/20
### Quantification of measured and modeled beampattern
* Beam center for individual click measurements `fig_single_bp_on_globe_20170920`:
	* if use _max point _after_ interpolation ('rx')_ as the beam center: Observed that the beam centers are either on top of the mic with max receiving level or between 2-3 mics with highest receicing levels. This means the beam center found in individual measured beampatterns are strongly affected by spatial under-sampling.
	* if use _center of the best-fitting ellipse ('ro')_ as the beam center, the locations of the beam center are more consistent with what one would expect by only looking at the -3 dB contours. Using this also potentially mitigates the bias from spatial under-sampling.
* Beam center for composite clicks: `fig_composite_click_avg_bp_20170920`
	* if use _max point _after_ interpolation ('rx')_ as the beam center: In the plot can see the beam center moves from outer side to inner side as frequency increases from 25 kHz to 45 kHz, but from 45 kHz to 50 kHz it's a _vertical downward_ movement. This pattern is reproduced for both left and right composite clicks.
	* if use _center of the best-fitting ellipse ('ro')_ as the beam center: the locations of the ellipse centers are generally fairly close to the max beam energy location
* Beam center for model results `fig_model_steer_h_bpctr`:
	* max beam energy location jumps around a bit across all steered angles, whereas center of best-fitting ellipse is very stable --> therefore use the latter in figures
* Tomorrow:
	* finish putting beam center on all figures
	* multi-freq beam center shifts for model (both individual and composite)
	* compare beam center location between data and model
	* compare el/az ratio between data and model
	* effects of missing one or more markers on the bat head



## 2017/09/21
### Cont: Quantification of measured and modeled beampattern
* Follow from yesterday: compare the max beam energy location vs averaged location for the top X%.
	* For composite measured clicks, the averaged location of all points >-1 dB normalized beam energy ('r^') is much closer to the location of the center of the best-fitting ellipse ('ro') than just the max beam energy location ('rx') (`fig_composite_click_avg_bp_20170920`)

	<img src=./img/fig_composite_click_avg_bp_20170920_batall_bin10_th0_35kHz_avg_bp_left.png width="400">	<img src=./img/fig_composite_click_avg_bp_20170920_batall_bin10_th0_35kHz_avg_bp_right.png width="400">
	* The above is the same for individual measured clicks: 'r^' is closer to 'ro' than 'rx' (`fig_single_bp_on_globe_20170920`)

	<img src=./img/fig_single_bp_on_globe_20170920_36134_02_c19_f35kHz_eckert4_rbf_mic.png width="400">	<img src=./img/fig_single_bp_on_globe_20170920_36134_02_c20_f35kHz_eckert4_rbf_mic.png width="400">
	* The above is also true for modeled beampattern: 'r^' is closer to 'ro' than 'rx' (`fig_model_steer_h_bpctr`)

	<img src=./img/fig_model_steer_h_bpctr_20160917_Ra-colony-rotear-0.5mm_x029t023_y000.0_z-05.0_2345_90deg_left_cntr.png width="800">	<img src=./img/fig_model_steer_h_bpctr_20160917_Ra-colony-rotear-0.5mm_x029t023_y000.0_z-05.0_2345_90deg_left_all.png width="800">

	* Based on the above results, the conclusion is that **it is safe to use the center of the best-fitting ellipse to denote the beam center**.
* Find beam center for composite modeled clicks:
	* It seems like the submitted paper used simulated results from 20161009, but should actually redo the simulation using rotation results from 20170308 (produced by `rotate_all_click_2010308`, note the missing `7` in the folder/file name). This would include the application of the newer version of `shift_rotate_bp` that does not kick out out-of-bnd points (which was corrected in `rotate_all_click_20161025`).
	* The new simulation code is `model_bp_proj_RaColony_20170921`, in which only 35 kHz projection was done, but parameters needed for projecting models at other frequencies based on beam center alignment at 35 kHz are saved. These critical parameters are `model.az_diff` and `model.el_diff`. Explanation of the saved variables is as follows:
		* `param`: parameters
		* `map`: definition of map system used
		* `data`: beam center calculated from data
		* `BP`: beam center calulated from beampattern model
		* `model_rot`: rotation results from projected beampattern
* **NOTE** All simulated beampattern in the original submission (v10 of fig folders) used mouth array position '3456'.
	* Considering using '2345' in the resubmission? The center of best-fitting ellipse and the averaged location for all points with normalized beam energy >-1 corresponds better to each other in '2345' scenario.
* **NOTE re. source locations used by original submission**
	*`model_bp_save_20161009` used a different numbering of sources that could cause confusion: used `SRC.idx_left([4,5,6,2]))` without the type of sorting used in `fig_model_steer_h_bpctr.m` now correct it to use the same type of sorting for source locations so that the code is more trackable --> but it turned out that model_bp_20161009 actually uses the source combination of '3456' correctly. Therefore results from `model_bp_proj_RaColony_20161025` used by the original submission were fine (indeed used sources '3456'.
	* `model_bp_save_multiple_20161014` used the same ordering of source locations without sorting like in `model_bp_save_20161009`, therefore results from `model_bp_proj_RaColony_multiple_20161025` used by the original submission were fine.



## 2017/09/22
### Cont: Quantification of measured and modeled beampattern
* Found out there were some problems with the code `model_bp_proj_RaColony_20170921` so spent quite some time to fix it. Now it seems fixed. Look for commit today (2017/09/22).
* Wrote a simple function `plot_bp_simple` to check interpolated beam pattern based on 1D or 2D data.
* Wrote `model_bp_proj_RaColony_diffonly_20170921` which only finds the critical parameters necessary for producing Monte Carlo simulation results, but leave the actual simulation to another code. The parameters obtained by this code are for enabling multi-freq Monte Carlo simulation. Here, _diffonly_ refers to the differences in beam center in the model and measurements. This is the most important parameter to get the simulated projection working.
* Wrote `model_bp_proj_RaColony_multifreq_20170921` which takes the output from `model_bp_proj_RaColony_diffonly_20170921` and simulate projected beam energy at mic locations at multiple frequencies. Need to revise `model_composite_RaColony3456_20170921` to work with output from this `...multifreq_20170921` code.



## 2017/09/23
### Cont: Quantification of measured and modeled beampattern
* Finished multi-freq model bp simulation and composite model bp compilation code:
  	          * `model_bp_proj_RaColony_diffonly_20170921`: calculate the shift between the beam center for measurements and model bp at 35 kHz.
		  * `model_bp_proj_RaColony_multifreq_20170921`: use the shift value from `...diffonly_20170921` to simulated model bp for 25:10:35 kHz.
		  * `model_composite_RaColony3456_20170921`: assemble simulated model bp values from `..._multifreq_20170921` to reconstructed composite model bp.
		  * the code above were executed in the order introduced.
* Wrote `add_bpctr_to_composite_output` to add beam center locations into composite beam measurement files. The routines are from `fig_composite_click_avg_bp_20170920`. This is to facilitate plotting the beam center in various representations (e.g., avg_bp or cntr).
* Compare the aspect ratio (el/az of best-fitting ellipse) between data and the phased array and piston model. Code is `fig_azel_distr_indiv_composite_model_20170923`. In here can see that the distributions of aspect ratio for the data and the phased array model are very similar:

	<img src=./img/fig_azel_distr_indiv_composite_model_20170923_batall_bin10_th0_nstd1.0_scatter.png width="600">    <img src=./img/fig_azel_distr_indiv_composite_model_20170923_batall_bin10_th0_nstd1.0_ar.png width="300">



## 2017/09/24
### Cont: Quantification of measured and modeled beampattern
* Need to re-run multi-freq model bp simulation code (`model_composite_RaColony3456_20170921.m`) because forgot to 60 kHz in the batch yesterday. The whole evaluated frequency range should be 20:5:60 kHz.
* Updated `model_composite_RaColony3456_20170921` to include getting -3dB contours and beam center for all frequencies. Different from the plotting routines for composite measured clicks, here all plots for model composite clicks are done in this code.
* Wrote `multifreq_indiv_bpctr_cntr_20170924` to plot the beam centers and -3dB contours at each frequency for all individual clicks (measured during exp). Note the beam center and contour locations were already aligned with the
* Revise `shift_rotate_bp.m`: adding upper limit for the time to try rotating and fitting ellipse (<2000) and varargout to track number of rotation (variable `rot_n`).
* Wrote a randomized model bp version of the multi-frequency Monte Carlo code: `model_bp_proj_RaColony_multiple_multifreq_20170924.m`. This code uses the same set of bp model as in `model_bp_proj_RaColony_multiple_20161025.m` but makes multi-frequency simulations.



## 2017/09/25
### Cont: Quantification of measured and modeled beampattern
* `check_head_marker_effect.m` checks the percentage and deviation of head aim angles between head aims derived from head markers vs from trajectories.
	* Analysis found 76% of the measurements were using head aim derived from head markers, and 24% from trajectories
	
	<img src=./img/check_head_marker_effect_goodidx.png width="1000">
	
	* Calculated the distribution of head aim differences (using clicks with which head markers are available to make this comparison), the result is as the following:
	
	<img src=./img/check_head_marker_effect_aim_diff.png width="400">
	
	* Make sure both codes below are working for extracting -3dB contours and beam center for all frequencies and for all clicks:
		* `multifreq_indiv_bpctr_cntr_20170924.m` --> data
		* `model_multifreq_indiv_bpctr_cntr_20170924.m` --> model
	* Comparing beam center location for data and model: mean
	
	<img src=./img/model_multifreq_indiv_bpctr_cntr_20170924_std1.0_azmean.png width="400">	<img src=./img/multifreq_indiv_bpctr_cntr_20170924_batall_azmean.png width="400">

	* Comparing beam center location for data and model: distribution **ectr**
	
	<img src=./img/model_multifreq_indiv_bpctr_cntr_20170924_std1.0_azdistr_ectr.png width="400">	<img src=./img/multifreq_indiv_bpctr_cntr_20170924_batall_azdistr_ectr.png width="400">

	* Comparing beam center location for data and model: distribution **top**
	
	<img src=./img/model_multifreq_indiv_bpctr_cntr_20170924_std1.0_azdistr_top.png width="400">	<img src=./img/multifreq_indiv_bpctr_cntr_20170924_batall_azdistr_top.png width="400">


 


* **NOTE** The `multifreq_composite_click_***_20170308` series of code starts by `multifreq_composite_click_20170308` and then calculate other related attributes, including `..._avg_bp` (beampattern for average click), `..._cntr` (multi-freq contour for average click), and `..._fit_elps` (fitting ellipse to -3dB contour at 35 kHz).




## TO-DO
* Use the plotting routines in `plot_indiv_click_rorate` to update those in the beampattern processing GUI, since doing conversion to map domain and plot using `contour` is much faster than calling `contourfm`



************************************************************************
## Codes used for paper figures

### Fig. 2-3
Experimental data, beampattern and multi-freq structure for both individual and composite clicks

- `fig_single_bp_on_globe_20160721`
- `fig_composite_click_avg_bp_20160808`
- `fig_el_cut_indiv_clicks_20160809` --> identical to `el_cut_indiv_clicks_20160722`, only file saving modifications
- `fig_multifreq_indiv_click_cntr_201600809`
- `fig_multifreq_composite_click_cntr_20160809`

Decided to use frequency spacing 25:5:55 kHz for multi-frequency plots. The reason is that the composite -3 dB contour for 20 kHz is very twisted for the threshold=0 and binsize=10deg case. The 20 kHz contour looks fine with higher threshold, but it's good to keep it consistent with the beampattern plot (Fig. 2A, 2B).

### Fig. 4

Code to plot mesh shapes:
- `fig_mesh_jigsaw` --> main code to plot mesh
- `draw_jigsaw_mesh` --> called by `fig_mesh_jigsaw` to draw mesh



************************************************************************
## List of BEM models calculated to date
Matching date/folder with model shape. Note all shapes are rotated such that the front faces +X axis

- **20160421**: bullethead (`bullethead_new_0.1.msh`), centered but not tilted
- **20160501**: bullethead (`bullethead_new_0.1.msh`), centered and tilted toward -Z-axis for -10 deg
- **20160727**: piston on square plate (`piston_new_0.1.msh`)
- **20160729**: circular piston shape (`piston_circ_0.05.msh`)
- **20160801**: bat head Ra\_224\_14k (`Rousettus_aegyptiacus_14k_0.88mm-remesh_fill_rot.msh`), centered and corrected for the head orientation
- **20160801b**: bat head Ra\_224\_7k (`Rousettus_aegyptiacus_7k_1.4mm-remesh_fill_rot.msh`), centered and corrected for the head orientation
- **20160803_Ra224-0.5mm**: Ra224 bat head, rotated around y-axis for -10 deg
- **20160803_Ra224-0.6mm**: Ra224 bat head, rotated around y-axis for -10 deg
- **20160811_Ra224-0.6mm**: Ra224 bat head, rotated around y-axis for -20 deg (mouth more parallel to x-axis)
- **20160812_Ra224-0.5mm**: Ra224 bat head, rotated around y-axis for -20 deg (mouth more parallel to x-axis)
- **20160812_bullethead-sc-0.5mm**: bullethead scaled to match dimension of Ra224, rotated around y-axis for 5 deg
- **20160814_bullethead-sc-0.6mm**: bullethead scaled to match dimension of Ra224, rotated around y-axis for 5 deg
- **20160814_Ra224-ear-0.5142mm**: Ra224 with ears, rotated around y-axis for -20 deg (mouth more parallel to x-axis)
- **20160817_Ra-colony-0.5mm**: Ra-colony bat head with ears (rotated to have mouth more parallel to x-axis)
- **20160817_Ra-colony-0.7mm**: Ra-colony bat head with ears (rotated to have mouth more parallel to x-axis)
- **20160817_bullethead-sc-colony-0.5mm**: bullethead scaled to match dimension of Ra-colony, rotated around y-axis for 5 deg
- **20160820_bullethead-sc-colony-0.5mm**: bullethead scaled to match dimension of Ra-colony, no rotation around y-axis
- **20160829_Ra-colony-0.7mm**: Ra-colony bat head with ears (rotated to have mouth more parallel to x-axis), this was calculated for many fine frequencies to verify if the anomaly is only at 25 kHz or also at nearby frequencies (--> the anomaly is only at 25 kHz)
- **20160829_Ra-colony-noear-0.5mm**: Ra-colony bat head *without* ears (rotated to have mouth more parallel to x-axis)
- <span style="background-color:yellow">**20160917_Ra-colony-rotear-0.5mm**</span>: Ra-colony bat head *with ears rotated slightly downward* (rotated to have mouth more parallel to x-axis) --> this slight ear rotation eliminated the 25 kHz anomaly --> this is the final mesh used in Fig. 4 in the paper

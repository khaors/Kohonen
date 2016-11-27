#########################################################################################################
#                                     Makefile for KOHONEN MAP
#########################################################################################################
#
# Comments
# - documentation requires robodoc and latex (pdflatex)
# - option -pg is to create info for gprof program
#
#########################################################################################################
FF = gfortran
FFLAGS=-c -g -pg -std=f2003 -Wall -fall-intrinsics -fbounds-check
OFLAGS=-g -pg -std=f2003 -Wall -fall-intrinsics -fbounds-check -o
LIBS= -llapack -lblas
RD = robodoc 
RDFLAGS1 = --src
RDFLAGS2 = --doc 
RDFLAGS3 = --singlefile --latex --sections --altlatex
PDF = pdflatex
########################################################################################################
#        Base classes
########################################################################################################

general_utilities.o: general_utilities.f90
	$(FF) $(FFLAGS) $<

kohonen_layer_parameters_utilities.o: kohonen_layer_parameters_utilities.f90
	$(FF) $(FFLAGS) $^

kohonen_prototype_utilities.o: kohonen_prototype_utilities.f90 distance_base_utilities.o
	$(FF) $(FFLAGS) $<

kohonen_pattern_utilities.o: kohonen_pattern_utilities.f90 kohonen_prototype_utilities.o
	$(FF) $(FFLAGS) $<

kohonen_map_base_utilities.o: kohonen_map_base_utilities.f90 kohonen_layer_parameters_utilities.o kohonen_pattern_utilities.o
	$(FF) $(FFLAGS) $<

self_organized_map_utilities.o: self_organized_map_utilities.f90 kohonen_layer_parameters_utilities.o kohonen_map_base_utilities.o kohonen_prototype_utilities.o distance_base_utilities.o factory_distance_utilities.o mtmod.o euclidean_distance_utilities.o manhattan_distance_utilities.o correlation_distance_utilities.o correlation2_distance_utilities.o itakura_saito_distance_utilities.o prediction_distance_utilities.o log_likelihood_distance_utilities.o wavelet_distance_utilities.o general_utilities.o 
	$(FF) $(FFLAGS) $<

distance_base_utilities.o: distance_base_utilities.f90
	$(FF) $(FFLAGS) $<

euclidean_distance_utilities.o: euclidean_distance_utilities.f90 distance_base_utilities.o
	$(FF) $(FFLAGS) $<

manhattan_distance_utilities.o: manhattan_distance_utilities.f90 distance_base_utilities.o
	$(FF) $(FFLAGS) $<

correlation_distance_utilities.o: correlation_distance_utilities.f90 distance_base_utilities.o
	$(FF) $(FFLAGS) $<

correlation2_distance_utilities.o: correlation2_distance_utilities.f90 distance_base_utilities.o
	$(FF) $(FFLAGS) $<

itakura_saito_distance_utilities.o: itakura_saito_distance_utilities.f90 distance_base_utilities.o
	$(FF) $(FFLAGS) $^

prediction_distance_utilities.o: prediction_distance_utilities.f90 distance_base_utilities.o
	$(FF) $(FFLAGS) $<

log_likelihood_distance_utilities.o: log_likelihood_distance_utilities.f90 distance_base_utilities.o
	$(FF) $(FFLAGS) $<

wavelet_distance_utilities.o: wavelet_distance_utilities.f90 distance_base_utilities.o
	$(FF) $(FFLAGS)  $<

dtw_distance_utilities.o: dtw_distance_utilities.f90 distance_base_utilities.o
	$(FF) $(FFLAGS)  $<

factory_distance_utilities.o: factory_distance_utilities.f90 euclidean_distance_utilities.o manhattan_distance_utilities.o correlation_distance_utilities.o correlation2_distance_utilities.o itakura_saito_distance_utilities.o prediction_distance_utilities.o log_likelihood_distance_utilities.o wavelet_distance_utilities.o dtw_distance_utilities.o distance_base_utilities.o
	$(FF) $(FFLAGS) $<

two_level_self_organized_map_utilities.o: two_level_self_organized_map_utilities.f90 kohonen_layer_parameters_utilities.o kohonen_map_base_utilities.o kohonen_prototype_utilities.o distance_base_utilities.o factory_distance_utilities.o mtmod.o
	$(FF) $(FFLAGS) $<

som_train_variables.o: som_train_variables.f90 kohonen_layer_parameters_utilities.o kohonen_map_base_utilities.o kohonen_prototype_utilities.o kohonen_pattern_utilities.o self_organized_map_utilities.o
	$(FF) $(FFLAGS) $<

som_predict_variables.o: som_predict_variables.f90 kohonen_layer_parameters_utilities.o kohonen_map_base_utilities.o kohonen_prototype_utilities.o kohonen_pattern_utilities.o self_organized_map_utilities.o
	$(FF) $(FFLAGS) $<

two_level_som_train_variables.o: two_level_som_train_variables.f90 kohonen_layer_parameters_utilities.o kohonen_map_base_utilities.o kohonen_prototype_utilities.o kohonen_pattern_utilities.o two_level_self_organized_map_utilities.o
	$(FF) $(FFLAGS) $^

mtmod.o: mtmod.f90
	$(FF) $(FFLAGS) $^

two_level_som_estimate_variables.o: two_level_som_estimate_variables.f90 kohonen_layer_parameters_utilities.o kohonen_map_base_utilities.o kohonen_prototype_utilities.o kohonen_pattern_utilities.o two_level_self_organized_map_utilities.o mtmod.o
	$(FF) $(FFLAGS) $^
########################################################################################################
#                                  Executable Files
########################################################################################################

all: som_train som_predict two_level_som_train two_level_som_estimate rsomlib

som_train: som_train.f90 som_train_variables.o self_organized_map_utilities.o kohonen_layer_parameters_utilities.o kohonen_map_base_utilities.o kohonen_prototype_utilities.o kohonen_pattern_utilities.o mtmod.o distance_base_utilities.o factory_distance_utilities.o euclidean_distance_utilities.o manhattan_distance_utilities.o correlation_distance_utilities.o correlation2_distance_utilities.o itakura_saito_distance_utilities.o prediction_distance_utilities.o log_likelihood_distance_utilities.o wavelet_distance_utilities.o dtw_distance_utilities.o general_utilities.o 
	$(FF) $^ $(LIBS) $(OFLAGS) $@

som_predict: som_predict.f90 som_predict_variables.o self_organized_map_utilities.o kohonen_layer_parameters_utilities.o kohonen_map_base_utilities.o kohonen_prototype_utilities.o kohonen_pattern_utilities.o distance_base_utilities.o factory_distance_utilities.o euclidean_distance_utilities.o manhattan_distance_utilities.o correlation_distance_utilities.o correlation2_distance_utilities.o itakura_saito_distance_utilities.o prediction_distance_utilities.o log_likelihood_distance_utilities.o wavelet_distance_utilities.o dtw_distance_utilities.o mtmod.o
	$(FF) $^  $(LIBS) $(OFLAGS) $@

two_level_som_train: two_level_som_train.f90 two_level_som_train_variables.o two_level_self_organized_map_utilities.o kohonen_layer_parameters_utilities.o kohonen_map_base_utilities.o kohonen_prototype_utilities.o kohonen_pattern_utilities.o distance_base_utilities.o factory_distance_utilities.o euclidean_distance_utilities.o manhattan_distance_utilities.o correlation_distance_utilities.o correlation2_distance_utilities.o itakura_saito_distance_utilities.o prediction_distance_utilities.o log_likelihood_distance_utilities.o wavelet_distance_utilities.o dtw_distance_utilities.o mtmod.o
	$(FF) $^ $(LIBS) $(OFLAGS) $@

two_level_som_estimate: two_level_som_estimate.f90 two_level_som_estimate_variables.o two_level_self_organized_map_utilities.o kohonen_layer_parameters_utilities.o kohonen_map_base_utilities.o kohonen_prototype_utilities.o kohonen_pattern_utilities.o distance_base_utilities.o factory_distance_utilities.o euclidean_distance_utilities.o manhattan_distance_utilities.o correlation_distance_utilities.o correlation2_distance_utilities.o itakura_saito_distance_utilities.o prediction_distance_utilities.o log_likelihood_distance_utilities.o wavelet_distance_utilities.o dtw_distance_utilities.o mtmod.o
	$(FF) $^ $(LIBS) $(OFLAGS) $@
	
rsomlib: self_organized_map_utilities.f90 kohonen_layer_parameters_utilities.f90 kohonen_map_base_utilities.f90 kohonen_prototype_utilities.f90 kohonen_pattern_utilities.f90 mtmod.f90 distance_base_utilities.f90 factory_distance_utilities.f90 euclidean_distance_utilities.f90 manhattan_distance_utilities.f90 correlation_distance_utilities.f90 correlation2_distance_utilities.f90 itakura_saito_distance_utilities.f90 prediction_distance_utilities.f90 log_likelihood_distance_utilities.f90 wavelet_distance_utilities.f90 dtw_distance_utilities.f90 general_utilities.f90
	rm *.o
	R CMD SHLIB self_organized_map_utilities.f90 kohonen_layer_parameters_utilities.f90 kohonen_map_base_utilities.f90 kohonen_prototype_utilities.f90 kohonen_pattern_utilities.f90 mtmod.f90 distance_base_utilities.f90 factory_distance_utilities.f90 euclidean_distance_utilities.f90 manhattan_distance_utilities.f90 correlation_distance_utilities.f90 correlation2_distance_utilities.f90 itakura_saito_distance_utilities.f90 prediction_distance_utilities.f90 log_likelihood_distance_utilities.f90 wavelet_distance_utilities.f90 dtw_distance_utilities.f90 general_utilities.f90 -fPIC -llapack -lblas
	rm *.o
###################################################################################################
####
###################################################################################################
clean:
	rm *.o
	rm *.mod
	rm *.pdf 
	rm *.tex
	rm kohonen_doc.*

###################################################################################################
####                                     Documentation
###################################################################################################
.PHONY: kohonen_doc
kohonen_doc:
	$(RD) --src /home/oscar/Oscar/MyProjects/Kohonen/ --doc $@ --singledoc --latex --sections --rc robodoc.rc --altlatex  
	$(PDF) $@.tex

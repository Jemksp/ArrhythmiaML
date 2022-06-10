out_dir = paste0("Directory containing all .txt tab delimited files for processing")
dir.create(out_dir, showWarnings = F)
setwd(out_dir)

############### Library##########################
library(caroline)
library(magrittr)
library(tidyr)
library(dplyr)
library(tibble)
library(openxlsx)
library(grid)
library(gridExtra)
library(ggplot2)
library(Matrix)
library(scatterplot3d)
library(rgl)
library(plotly)
library(gridExtra)
library(gtable)
library(lattice)
library(stringr)

############### Params###############
#Primary Parameters (Sensitivity ~ 0.05 to 0.3, lsfactor ~ 1.1-1.35, rsfactor = 1, start = at least 4)
sensitivity = 0.8
lsfactor = 1
rsfactor = 1
start = 4

#threshold for differential detection (0.01 - 0.5) (0.2)
dethreshold = 0.4
#threshold to drop peaks if no minima is in this search space (0.5 - 2) (0.8)
minsearch = 2
#threshold to drop minima if more than 1 in this search space (0 - 2) (1)
minsearch2 = 0.5
#threshold to drop minima if no peaks is in this search space (*1-3) (Should be the same as minsearch)
minsearch3 = minsearch * 1

#thresholds for bleaching (Bleaching = 0 for no bleaching, 1 for bleaching present) 
bleaching = 0
#threshold to remove peaks below gating line (0 to turn off, else between 0.25-0.7)
gstrict = 0.25
#threshold for gate considering bleaching (More than 1-fold drop in intensity) (0.1 ~ 0.4)
gstrict2 = 0.4
#threshold for gate to detect bleaching (0.2 ~ 0.7) (A larger number triggers bleaching detection more easily)
gstrict3 = 0.6
#If (1), Set threshold to (Maximum Intensity + Min intensity) * gstrict 
gstrict4 = 2
#Column classifier (Healthy = 0, Arrhythmic = 1)
Arrhythmia_classifier = 0

#Debugging (Run it on 0 to remove all the troubleshooting lines)
testbug = 1
gatetail = 6 #Leave on 6, number of peaks at the end to consider for bleaching

#################### Run####

for (hide in 1) {
  
  #Input file to be analysed ----
  defined_reads = list.files()
  
  Classifier = "Arrhythmia_status"
  
  
  #Preprocess ----
  tencol <- function(vector_input) {
    if (length(vector_input) %% 10 != 0) {
      vector_input <- c(vector_input, vector(mode = "character", length = (10 - length(vector_input) %% 10)))
    } else {
      vector_input
    }
  }
  
  overall_output = data.frame()
  overall_graph = list()
  raster_count = 0
  overall_raster = data.frame()
  graph_count = 1
  ROI_graph = list()
  
  remove = vector()
  for (i in (1:length(defined_reads))) {
    if (substring(defined_reads[i], nchar(defined_reads[i])-2,nchar(defined_reads[i])) == "txt") {
      remove = c(remove,i)
    }
  }
  defined_reads = defined_reads[remove]
  
  #Iteration for all samples and ROIs ----
  for (sample in 1:length(defined_reads)) {
    
    raw_data = read.delim(file = defined_reads[sample], fileEncoding = "UTF-16LE")
    df <- raw_data[-c(3:5)]
    colnames(df) <- c("Name","Time", matrix("ROI ", ncol = ncol(df)-2) %>%
                        paste0(1:(ncol(df)-2)))
    pdf_outputlist = list()
    outputh = vector()
    
    #Compiling roi
    for (roi in 3:(ncol(df))) {
    #for (roi in 5) { #troubleshooting
      page_title = as.character(df$Name[1]) %>%
        str_remove(".nd2")
      roi_title = colnames(df)[roi]
      iegate = 0
      
      raster_count = raster_count + 1
      dfc = data.frame(df$Time,df[roi])
      colnames(dfc) <- c("Time", "Intensity")
      Intensity = df[[roi]]
      
      #Sensitivity optimisation
      if (substring(defined_reads[sample],1,3) != "ORG") {
        sense1 = df$Time - df$Time[1]
        sense = floor(min(which(sense1>sensitivity)))
        ls = floor((lsfactor-1)*sense)
        rs = floor((rsfactor-1)*sense)
        
        #Computing differential
        de_df = vector()
        for (i in 1:(length(df$Time)-1)) {
          de_df[i] = dfc$Intensity[i+1] - dfc$Intensity[i]
        }
        
        de_t = mean(de_df[order(-de_df)][1:10]) * dethreshold
        de_outs = vector()
        for (i in start:(length(de_df)-1)) {
          if ((de_df[i-2] <= de_t | de_df[i-3] <= de_t) &
              de_df[i-1] <= de_t &
              de_df[i] >= de_t) {
            de_outs = c(de_outs, i)
          }
        }
        
        #Dropping doubled minimas
        de_debug_1 = de_outs
        de_t2 = vector()
        de_t2_count = 1
        for (i in 2:length(de_outs)) {
          de_double = de_outs - de_outs[i-1]
          if (any(between(de_double,floor(-sense*minsearch2),-1))) {
            de_t2_test1 = which(between(de_double,floor(-sense*minsearch2),-1))
            if (length(de_t2_test1 > 1)) {
              de_t2_test1 = de_t2_test1[length(de_t2_test1)]
            }
            de_t2_test2 = de_t2_test1 +1
            de_t2_test3 = dfc$Intensity[de_outs[de_t2_test1]] - dfc$Intensity[de_outs[de_t2_test2]]
            # print(de_t2_test1)
            if (de_t2_test3 < 0) {
              de_t2[de_t2_count] = de_t2_test2
            } 
            else {
              de_t2[de_t2_count] = de_t2_test1
            }
            de_t2_count = de_t2_count + 1
          }
        }
        if (de_t2_count != 1) {
        de_outs = de_outs[-de_t2]
        }
        
        
        # #Testing
        # de_outs_t = de_outs
        # de_double = vector()
        # for (j in 2:length(de_outs_t)) {
        #   de_double[j-1] = de_outs_t[j] - de_outs_t[j-1]
        # }
        # remove = which(de_double < floor(sense*minsearch2))
        # de_double2 = vector()
        # for (k in 2:length(remove)) {
        #   de_double2[k-1] = remove[k] - remove[k-1]
        # }
        #   #remove[length(de_double2) + 1]
        
        
        #Computing Peaks
        outs = vector()
        for (i in start:(length(df$Time)-(sense+rs))) {
          #if intensity point is the maximum of sense values on the left and right & not equal to another point of sense value to the left, record its value
          if (Intensity[i] == max(Intensity[(if((i-sense-ls) < 1) {1} else {i-sense-ls}):(if((i+sense+rs) > length(df$Time)) {length(df$Time)} else {i+sense+rs})]) & Intensity[i] != max(Intensity[(if((i-sense) < 1) {1} else {i-sense}) :(i-1)]))   {
            outs = c(outs, i)
          }
        }
        
        #Computing troughs
        touts = vector()
        for (i in 1:(length(outs) - 1)) {
          touts = c(touts, (outs[i]:outs[i+1])[which.min(Intensity[outs[i]:outs[i+1]])])
        }
        touts_debug_original = touts
        
        #Selecting higher peaks
        outs_debug_1 = outs
        gate = max(dfc$Intensity[tail(outs,gatetail)]) - min(dfc$Intensity[touts[1:2]])
        gate3 = gstrict3*(max(dfc$Intensity[outs[1:3]]) - min(dfc$Intensity[touts[1:2]]))
        outs_debug_original <- outs
        if ((gate > gate3 & gstrict > 0) | gstrict4 == 2) {
          threshold = min(dfc$Intensity[touts[1:2]]) + gstrict*gate
          if (gstrict4 > 0) {
            threshold = min(dfc$Intensity) + ((max(dfc$Intensity) - min(dfc$Intensity)) * gstrict)
          }
          gated = which(dfc$Intensity[outs] > threshold)
          outs2 <- outs[gated]
          
          if (bleaching == 1) {
          midpoint <- floor(length(outs2)/4)
          midpointint <- max(dfc$Intensity[outs2[midpoint:(midpoint+5)]])
          gate2 = midpointint - min(dfc$Intensity[touts[1:2]])
          threshold2 = min(dfc$Intensity[touts[1:2]]) + gstrict*gate2
          gated2 = c(which(dfc$Intensity[outs2[1:(midpoint+5)]] > threshold2), which(outs2 > outs2[midpoint+5]))
          outs <- outs2[gated2]
          }
          
          if (bleaching == 0) {
            threshold2 = threshold
            outs <- outs2
          }
          
          iegate = 1
        }
        else {if (gstrict >0)
        {
          #Things to do: make gstrict a gradient as well
          gradient_drop = min(dfc$Intensity[1:outs[2]]) - min(dfc$Intensity[(length(dfc$Time)/2):length(dfc$Time)])
          gradient_depolarization = max(dfc$Intensity[outs]) - min(dfc$Intensity[touts[1:2]])
          gradient_start_threshold = min(dfc$Intensity[1:outs[2]]) + gstrict2*gradient_depolarization
          gradient_drop_unit = gradient_drop / tail(dfc$Time,1)
          outs2 = vector()
          gradient_count = 1
          for (i in 1:length(outs)) {
            outs2[i] = dfc$Intensity[outs][i] - (gradient_start_threshold - (dfc$Time[outs][i] * gradient_drop_unit))
          }
          gradient_test = which(outs2 < 0)
          if (length(gradient_test > 0)) {
            outs <- outs[-gradient_test]
          }
          gradient_plot = data.frame(Intensity = c(gradient_start_threshold, gradient_start_threshold - (tail(dfc$Time,1)*gradient_drop_unit)),
                                     Time = c(0, tail(dfc$Time,1)))
          iegate = 2
        }
        }
        
        #Recalibrate sense
        sense2 = sense
        sense_outs = vector()
        for (i in 1:(length(outs) - 1)) {
          sense_outs[i] = outs[i+1] - outs[i]
        }
        sense_ave = mean(sense_outs)
        sense = floor(sense_ave)
        
        #Selecting peaks with detected minima
        outs_debug_2 = outs
        outs3 = vector()
        de_out_count = 1
        for (i in 1:length(outs)) {
          de_test = de_outs - outs[i]
          if (any(between(de_test, -floor(sense*minsearch), -1))) {
            outs3[de_out_count] = outs[i]
            de_out_count = de_out_count+1
          }
        }
        outs <- outs3
        
        #Selecting minima with detected peaks (Every peak can only keep 1 minima)
        de_debug_2 = de_outs
        de_outs2 = vector()
        de_out_count = 1
        de_outs3 = de_outs
        for (i in 1:length(outs)) {
          de_test = de_outs3 - outs[i]
          if (any(between(de_test, -floor(sense*minsearch3), -1))) {
            selection = which(between(de_test, -floor(sense*minsearch3), -1))
            if (length(selection) > 1) {
              de_outs3[1:(tail(selection,1))] = -999
              selection_index = which.min(dfc$Intensity[de_outs[selection]])
              selection = selection[selection_index]
            } 
            if (length(selection == 1)) {
              de_outs3[1:selection] = -999
            }
            de_outs2[de_out_count] = selection
            de_out_count = de_out_count + 1
          }
        }
        de_outs <- de_outs[de_outs2]
        
        #Pairing minimas to peaks 1 to 1 (Every minima can only keep 1 peak, drop other peaks)
        outs_debug_3 = outs
        outs4 = vector()
        de_out_count = 1
        for (i in 1:length(de_outs)) {
          de_test = outs - de_outs[i]
          outs4[de_out_count] = head(which(de_test > 0),1)
          de_out_count = de_out_count + 1
        }
        outs = outs[outs4]

        #Plotting minima
        df_min <- df %>%
          select (2, all_of(roi)) %>%
          slice(de_outs)
        colnames(df_min)[2] <- "minima"
        
        #Finalizing peaks
        peaks <- df$Time[outs]
        
        # Troubleshooting
        # de_debug_1 #original
        # de_debug_2 #after doubled minima removal
        # de_outs #after no peak removal
        # outs_debug_1 #original
        # outs_debug_2 #after gstrict
        # outs_debug_3 #after no minima removal
        # outs #after no 1:1 minima removal
        
        #Recompute Troughs
        touts = vector()
        for (i in 1:(length(outs) - 1)) {
          touts = c(touts, (outs[i]:outs[i+1])[which.min(Intensity[outs[i]:outs[i+1]])])
        }

        
        #For plotting peaks on graphs
        dat <- df %>% 
          select (2, all_of(roi)) %>% 
          slice(outs)
        colnames(dat)[2] <- "peakpoint"
      }
      
      #Computing wavelengths
      wavelength = 1:((length(peaks)-1))
      for (i in 1:(length(peaks)-1)) {
        wavelength[i] = peaks[i+1] - peaks[i]
      }
      
      #Compute Peak and Trough fall intensities and threshold for APD
      #Take the Peak Intensity minus the Minima intensity of the next peak
      apd100_int <- dfc$Intensity[outs[-length(outs)]] - dfc$Intensity[touts]
      apd30_int <- apd100_int*0.3
      apd60_int <- apd100_int*0.6
      apd90_int <- apd100_int*0.9
      apd30 = vector()
      apd60 = vector()
      apd90 = vector()
      apd30v = vector()
      apd60v = vector()
      apd90v = vector()
      apd30b60 = vector()
      apd30b90 = vector()
      apd60b90 = vector()
      cad30 = vector()
      cad60 = vector()
      cad90 = vector()
      cad30b60 = vector()
      cad30b90 = vector()
      cad60b90 = vector()
      # cad30v = vector()
      # cad60v = vector()
      # cad90v = vector()
      #repeats = length(head(sort(c(outs,touts)),-1))/2
      repeats = length(apd100_int)
      peak_integral = vector()
      #Compute APDurations, instability and triangulation mean and SD
      for (i in 1:repeats) {
        #Compute all intensities from a peak to next minima, which intensity value is the first point where it drops below 30%,60%,90%, Find the time value and subtract the time of the peak
        # group_pt = Intensity[head(sort(c(outs,touts)),-1)[(i*2)-1]:head(sort(c(outs,touts)),-1)[i*2]]
        group_pt = dfc$Intensity[outs[i]:touts[i]]
        group_pt2 = group_pt[1] - group_pt
        #Repolarization Vector index values - Which df index goes below the repolarization % mark of that peak
        apd30v[i] = outs[i] + head(which(group_pt2 > apd30_int[i]),1) - 1
        apd60v[i] = outs[i] + head(which(group_pt2 > apd60_int[i]),1) - 1
        apd90v[i] = outs[i] + head(which(group_pt2 > apd90_int[i]),1) - 1
        
        #Time values
        apd30[i] = dfc$Time[apd30v[i]] - dfc$Time[outs[i]]
        apd60[i] = dfc$Time[apd60v[i]] - dfc$Time[outs[i]]
        apd90[i] = dfc$Time[apd90v[i]] - dfc$Time[outs[i]]
        
        #Fractions
        apd30b60[i] = apd30[i]/apd60[i]
        apd30b90[i] = apd30[i]/apd90[i]
        apd60b90[i] = apd60[i]/apd90[i]
        
        #Calcium Transient Durations
        cad30[i] = dfc$Time[apd30v[i]] - dfc$Time[de_outs[i]]
        cad60[i] = dfc$Time[apd60v[i]] - dfc$Time[de_outs[i]]
        cad90[i] = dfc$Time[apd90v[i]] - dfc$Time[de_outs[i]]
        
        #Fractions
        cad30b60[i] = cad30[i]/cad60[i]
        cad30b90[i] = cad30[i]/cad90[i]
        cad60b90[i] = cad60[i]/cad90[i]
        
        #Computing Peak Integrals
        if (length(outs) > 2) {
            frames_for_integral = de_outs[i]:apd90v[i]
            integral_sum = 0
            for (j in 1:length(frames_for_integral)) {
              integral_sum = integral_sum + 
                (
                  (dfc$Intensity[j+1] - 
                     (df$Time[frames_for_integral[j]] - df$Time[frames_for_integral[1]]) *
                     ((dfc$Intensity[tail(frames_for_integral,1)] - dfc$Intensity[frames_for_integral[1]])/(df$Time[tail(frames_for_integral,1)] - df$Time[frames_for_integral[1]])) +
                     (Intensity[frames_for_integral[1]])                                    
                  )
                  *(df$Time[j+1] - df$Time[j])
                )
              #Integration of the peak using 
              #Summation of (intensity[point+1] - pseudo_base_intensity)
              ##calculated using y = mx+c whereby m = gradient from Int[1] oto Int[APD90], x is Time distance of point from [1], and C is Int[base]
              #multiplied with the change in time between point+1 and point
            }
            peak_integral[i] = integral_sum
        }
        
        # apd30[i] = dfc$Time[min(which(group_pt < (group_pt[1] - apd30_int[i])))+(head(sort(c(outs,touts)),-1)[(i*2)-1] - 1)] - dfc$Time[head(sort(c(outs,touts)),-1)[(i*2)-1]]
        # apd60[i] = dfc$Time[min(which(group_pt < (group_pt[1] - apd60_int[i])))+(head(sort(c(outs,touts)),-1)[(i*2)-1] - 1)] - dfc$Time[head(sort(c(outs,touts)),-1)[(i*2)-1]]
        # apd90[i] = dfc$Time[min(which(group_pt < (group_pt[1] - apd90_int[i])))+(head(sort(c(outs,touts)),-1)[(i*2)-1] - 1)] - dfc$Time[head(sort(c(outs,touts)),-1)[(i*2)-1]]
        # #Vector index values
        # apd30v[i] = min(which(group_pt < (group_pt[1] - apd30_int[i])))+(head(sort(c(outs,touts)),-1)[(i*2)-1] - 1)
        # apd60v[i] = min(which(group_pt < (group_pt[1] - apd60_int[i])))+(head(sort(c(outs,touts)),-1)[(i*2)-1] - 1)
        # apd90v[i] = min(which(group_pt < (group_pt[1] - apd90_int[i])))+(head(sort(c(outs,touts)),-1)[(i*2)-1] - 1)
      }
      # for (i in 2:repeats) {
      #   group_pt = Intensity[head(sort(c(outs,touts)),-1)[(i*2)-1]:head(sort(c(outs,touts)),-1)[i*2]]
      #   #CAD Time values
      #   cad30[i-1] = dfc$Time[min(which(group_pt < (group_pt[1] - apd30_int[i])))+(head(sort(c(outs,touts)),-1)[(i*2)-1] - 1)] - dfc$Time[touts[i-1]]
      #   cad60[i-1] = dfc$Time[min(which(group_pt < (group_pt[1] - apd60_int[i])))+(head(sort(c(outs,touts)),-1)[(i*2)-1] - 1)] - dfc$Time[touts[i-1]]
      #   cad90[i-1] = dfc$Time[min(which(group_pt < (group_pt[1] - apd90_int[i])))+(head(sort(c(outs,touts)),-1)[(i*2)-1] - 1)] - dfc$Time[touts[i-1]]
      #   #CAD Vector indexes
      #   # cad30v[i-1] = min(which(group_pt < (group_pt[1] - apd30_int[i])))+(head(sort(c(outs,touts)),-1)[(i*2)-1] - 1)
      #   # cad60v[i-1] = min(which(group_pt < (group_pt[1] - apd60_int[i])))+(head(sort(c(outs,touts)),-1)[(i*2)-1] - 1)
      #   # cad90v[i-1] = min(which(group_pt < (group_pt[1] - apd90_int[i])))+(head(sort(c(outs,touts)),-1)[(i*2)-1] - 1)
      #   #Fractions
      #   cad30b60[i-1] = cad30[i-1]/cad60[i-1]
      #   cad30b90[i-1] = cad30[i-1]/cad90[i-1]
      #   cad60b90[i-1] = cad60[i-1]/cad90[i-1]
      # }
      
      #Compilations
      df_apd = data.frame(head(peaks,-1), apd30, apd60, apd90, apd30b60, apd30b90, apd60b90)
      colnames(df_apd) = c("Peak", "APD30", "APD60", "APD90","APD30by60","APD30by90","APD60by90")
      df_cad = data.frame(df_apd[[1]], cad30, cad60, cad90, cad30b60, cad30b90, cad60b90)
      colnames(df_cad) = c("Peak", "CAD30", "CAD60", "CAD90","CAD30by60","CAD30by90","CAD60by90")
      
      apd30_sd = sd(df_apd[[2]])
      apd60_sd = sd(df_apd[[3]])
      apd90_sd = sd(df_apd[[4]])
      
      #Graphing APD30,60,90
      apdg <- df %>% 
        select (2, all_of(roi)) %>% 
        slice(c(apd30v,apd60v,apd90v))
      colnames(apdg)[2] <- "APD"
      
      #Accounting for bleaching to measure intensity variations
      processed_peak_intensity = head(dfc$Intensity[outs] - dfc$Intensity[de_outs],-1)
      peak_intensity_sd = sd(processed_peak_intensity)
      
      #Computing peak integral variables
      # peak_integral = vector()
      # if (length(touts) > 1) {
      #   for (i in 1:(length(touts)-1)) {
      #     frames_for_integral = touts[i]:apd90v[i+1]
      #     integral_sum = 0
      #     for (j in 1:length(frames_for_integral)) {
      #       integral_sum = integral_sum + 
      #         (
      #         (Intensity[j+1] - 
      #            (df$Time[frames_for_integral[j]] - df$Time[frames_for_integral[1]]) *
      #            ((Intensity[tail(frames_for_integral,1)] - Intensity[frames_for_integral[1]])/(df$Time[tail(frames_for_integral,1)] - df$Time[frames_for_integral[1]])) +
      #            (Intensity[frames_for_integral[1]])                                    
      #         )
      #         *(df$Time[j+1] - df$Time[j])
      #         )
      #       #Integration of the peak using 
      #       #Summation of (intensity[point+1] - pseudo_base_intensity)
      #       ##calculated using y = mx+c whereby m = gradient from Int[1] oto Int[APD90], x is Time distance of point from [1], and C is Int[base]
      #       #multiplied with the change in time between point+1 and point
      #     }
      #     peak_integral[i] = integral_sum
      #   }
      # }
      peak_integral_mean = mean(peak_integral)
      peak_integralnm_mean = mean(peak_integral/processed_peak_intensity)
      peak_integralfr_mean = peak_integral_mean/mean(wavelength)
      peak_integralnmfr_mean = peak_integralnm_mean/mean(wavelength)
      peak_integral_sd = sd(peak_integral)
      peak_integralnm_sd = sd(peak_integral/processed_peak_intensity)
      peak_integralb90_mean = mean(peak_integral/cad90)
      peak_integralb90_sd = sd(peak_integral/cad90)
      
      #Computing Rises
      crise = df$Time[outs] - df$Time[de_outs]
      
      #Rises by CAD90
      riseb90 <- crise[-length(crise)] / cad90
      
      #Rises by Intensity (Depolarization Intensity)
      risevelocity <- (dfc$Intensity[outs] - dfc$Intensity[de_outs])/crise
      
      #Compiling peaks into 10 values per row
      #df_peaks = t(matrix(tencol(peaks), nrow = 10))
      
      #Compiling wavelengths
      # df_wavelength = t(matrix(tencol(wavelength), nrow = 10))
      # df_wlmean = round(mean(wavelength),4)
      # df_wlsdev = round(sd(wavelength),4)
      
      #Compiling Rises
      #df_rises = t(matrix(tencol(crise), nrow = 10))
      # df_meanrise = round(mean(crise),4)
      # df_sdevrise = round(sd(crise),4)
      
      #Computing APDfr/frsd, Triangulationfr/frsd (fraction of mean wavelength)
      apd30fr = apd30/(mean(wavelength))
      apd30fr_mean = mean(apd30fr)
      apd30fr_sd = sd(apd30fr)
      
      apd60fr = apd60/(mean(wavelength))
      apd60fr_mean = mean(apd60fr)
      apd60fr_sd = sd(apd60fr)
      
      apd90fr = apd90/(mean(wavelength))
      apd90fr_mean = mean(apd90fr)
      apd90fr_sd = sd(apd90fr)
      
      triangulation_mean = mean(df_apd[[4]] - df_apd[[2]])
      triangulation_sd = sd(df_apd[[4]] - df_apd[[2]])
      triangulation_frmean = mean((df_apd[[4]] - df_apd[[2]])/(mean(wavelength)))
      triangulation_frsd = sd((df_apd[[4]] - df_apd[[2]])/(mean(wavelength)))
      
      cad30fr = cad30/(mean(wavelength))
      cad30fr_mean = mean(cad30fr)
      cad30fr_sd = sd(cad30fr)
      
      cad60fr = cad60/(mean(wavelength))
      cad60fr_mean = mean(cad60fr)
      cad60fr_sd = sd(cad60fr)
      
      cad90fr = cad90/(mean(wavelength))
      cad90fr_mean = mean(cad90fr)
      cad90fr_sd = sd(cad90fr)
      
      #Overall spreadsheet output
      single_output = data.frame(Arrhythmia_classifier, index = 1,
                                 mean(wavelength), sd(wavelength), mean(crise), sd(crise), mean(riseb90), sd(riseb90), 
                                 mean(df_apd[[2]]), mean(df_apd[[3]]), mean(df_apd[[4]]), apd30_sd, apd60_sd, apd90_sd,
                                 mean(df_apd[[5]]), mean(df_apd[[6]]), mean(df_apd[[7]]), sd(df_apd[[5]]), sd(df_apd[[6]]), sd(df_apd[[7]]), 
                                 apd30fr_mean, apd60fr_mean, apd90fr_mean, apd30fr_sd, apd60fr_sd, apd90fr_sd, 
                                 mean(df_cad[[4]]), sd(df_cad[[4]]), cad90fr_mean, cad90fr_sd, mean(risevelocity), sd(risevelocity),
                                 peak_integral_mean, peak_integral_sd, peak_integralb90_mean, peak_integralb90_sd, peak_integralfr_mean, 
                                 cad30fr_mean, cad60fr_mean, cad30fr_sd, cad60fr_sd, 
                                 mean(df_cad[[2]]), mean(df_cad[[3]]), sd(df_cad[[2]]), sd(df_cad[[3]]), 
                                 mean(df_cad[[5]]), mean(df_cad[[6]]), mean(df_cad[[7]]),  sd(df_cad[[5]]), sd(df_cad[[6]]), sd(df_cad[[7]]),
                                 peak_integralnm_mean, peak_integralnmfr_mean, peak_integralnm_sd,  peak_intensity_sd, 
                                 triangulation_mean, triangulation_sd, triangulation_frmean, triangulation_frsd, 
                                 row.names = paste(page_title, roi_title))
      overall_output = rbind(overall_output, single_output)
      
      single_raster = data.frame(peaks, raster_count, sample, roi-2, page_title, c(wavelength,""))
      #Use this to create the required stuff for poincare plots
      overall_raster = rbind(overall_raster, single_raster)
      
      #Arranging grobs in 1 page ------------------
      # g1 = tableGrob(df_peaks)
      # title <- textGrob("Beats",gp=gpar(fontsize=25))
      # padding <- unit(5,"mm")
      # tg1 <- gtable_add_rows(g1, heights = grobHeight(title) + padding, pos = 0)
      # tg1 <- gtable_add_grob(tg1, title, 1, 1, 1, ncol(tg1))
      # 
      # g2 = tableGrob(df_wavelength)
      # title <- textGrob("Beat to Beat duration",gp=gpar(fontsize=25))
      # padding <- unit(5,"mm")
      # tg2 <- gtable_add_rows(g2, heights = grobHeight(title) + padding, pos = 0)
      # tg2 <- gtable_add_grob(tg2, title, 1, 1, 1, ncol(tg2))
      # 
      # g3 = tableGrob(df_meanrise, theme = ttheme_default(base_size = 24))
      # title <- textGrob("Depolarization duration Mean",gp=gpar(fontsize=21))
      # padding <- unit(5,"mm")
      # tg3 <- gtable_add_rows(g3, heights = grobHeight(title) + padding, pos = 0)
      # tg3 <- gtable_add_grob(tg3, title, 1, 1, 1, ncol(tg3), clip = 'off')
      # 
      # g4 = tableGrob(df_sdevrise, theme = ttheme_default(base_size = 24))
      # title <- textGrob("Depolarization duration SD",gp=gpar(fontsize=21))
      # padding <- unit(5,"mm")
      # tg4 <- gtable_add_rows(g4, heights = grobHeight(title) + padding, pos = 0)
      # tg4 <- gtable_add_grob(tg4, title, 1, 1, 1, ncol(tg4), clip = 'off')
      # 
      # g5 = tableGrob(df_wlmean, theme = ttheme_default(base_size = 24))
      # title <- textGrob("Beat to Beat Mean",gp=gpar(fontsize=21))
      # padding <- unit(5,"mm")
      # tg5 <- gtable_add_rows(g5, heights = grobHeight(title) + padding, pos = 0)
      # tg5 <- gtable_add_grob(tg5, title, 1, 1, 1, ncol(tg5), clip = 'off')
      # 
      # g6 = tableGrob(df_wlsdev, theme = ttheme_default(base_size = 24))
      # title <- textGrob("Beat to Beat SD",gp=gpar(fontsize=21))
      # padding <- unit(5,"mm")
      # tg6 <- gtable_add_rows(g6, heights = grobHeight(title) + padding, pos = 0)
      # tg6 <- gtable_add_grob(tg6, title, 1, 1, 1, ncol(tg6), clip = 'off')
      
      if (testbug == 1) {
      if (iegate == 1) {
        gplot = ggplot() + 
        geom_line(data = dfc, aes(x = Time, y = Intensity)) +
        #geom_point(data = dfc, aes(x = Time, y = Intensity), color = 'gray', size = 0.1) + 
        geom_point(data = dat, aes(x = Time, y = peakpoint), color = 'red', size = 2.5) +
        geom_point(data = apdg, aes(x = Time, y = APD), color = 'blue', size = 0.9) +
        geom_point(data = df_min, aes(x = Time, y = minima), color = 'orange', size = 2) +
        geom_hline(yintercept=threshold, color = 'springgreen3') +
        geom_hline(yintercept=threshold2, color = 'springgreen3') +
        labs(title = paste(page_title, "ROI =", roi-2 ,"SLRS,DMMM,BGG =", sensitivity,lsfactor,rsfactor,start,",",dethreshold,minsearch,minsearch2,minsearch3,",",bleaching,gstrict,gstrict2)) +
        theme_bw() +
        xlab("Time [S]") +
        ylab("Intensity")
        
        if (bleaching == 1) {
          gplot = gplot + geom_vline(xintercept=dfc$Time[outs2[midpoint]], color = 'springgreen3')
        }
      }
      if (iegate == 2) {
        gplot = ggplot() + 
          geom_line(data = dfc, aes(x = Time, y = Intensity)) +
          #geom_point(data = dfc, aes(x = Time, y = Intensity), color = 'gray', size = 0.1) + 
          geom_point(data = dat, aes(x = Time, y = peakpoint), color = 'red', size = 2.5) +
          geom_point(data = apdg, aes(x = Time, y = APD), color = 'blue', size = 0.9) +
          geom_point(data = df_min, aes(x = Time, y = minima), color = 'orange', size = 2) +
          geom_line(data = gradient_plot, aes(x = Time, y = Intensity), color = 'springgreen3') +
          labs(title = paste(page_title, "ROI =", roi-2 ,"SLRS,DMMM,BGG =", sensitivity,lsfactor,rsfactor,start,",",dethreshold,minsearch,minsearch2,minsearch3,",",bleaching,gstrict,gstrict2)) +
          theme_bw() +
          xlab("Time [S]") +
          ylab("Intensity")
      }
      }
      if (testbug == 0 | gstrict == 0) {
      gplot = ggplot() +
        geom_line(data = dfc, aes(x = Time, y = Intensity)) +
        #geom_point(data = dfc, aes(x = Time, y = Intensity), color = 'gray', size = 0.1) +
        geom_point(data = dat, aes(x = Time, y = peakpoint), color = 'red', size = 2.5) +
        geom_point(data = apdg, aes(x = Time, y = APD), color = 'blue', size = 0.9) +
        geom_point(data = df_min, aes(x = Time, y = minima), color = 'orange', size = 1.5) +
        labs(title = paste(page_title, "ROI =", roi-2 ,"SLRS,DMMM,BGG =", sensitivity,lsfactor,rsfactor,start,",",dethreshold,minsearch,minsearch2,minsearch3,",",bleaching,gstrict,gstrict2,gstrict3)) +
        theme_bw() +
        xlab("Time [S]") +
        ylab("Intensity")
      }
      
      # lay <- rbind(c(1,1,1,1,3,3),
      #              c(1,1,1,1,4,4),
      #              c(2,2,2,2,5,5),
      #              c(2,2,2,2,6,6),
      #              c(7,7,7,7,7,7),
      #              c(7,7,7,7,7,7),
      #              c(7,7,7,7,7,7),
      #              c(7,7,7,7,7,7))
      
      # gout = grid.arrange(tg1, tg2, tg3, tg4, tg5, tg6, gplot, layout_matrix = lay, top = textGrob(paste(page_title,roi_title, "S =", sensitivity), gp = gpar(fontsize = 18)))
      # 
      # #Compiling output lists
      # pdf_outputlist[[roi - 2]] <- gout
      ROI_graph[[graph_count]] <- gplot
      #end of 1 roi
      graph_count = graph_count + 1
    }
    
    #Output individual PDF ----------------------------
    # ggsave(str_replace(defined_reads[sample], ".txt", ".pdf"), 
    #        marrangeGrob(pdf_outputlist, nrow = 1, ncol = 1, top = textGrob(" ", gp = gpar(fontsize = 10))), 
    #        width = 15, height = 14, limitsize = F)
    # dev.off()

    #end of 1 sample
  }
  
  
  #Output compilations -----------------------
  #Raster compilation (Poincare plot)
  colnames(overall_raster) <- c("Peaks", "Total Count", "Sample Count", "ROI Count", "Doc Title", "RR_Interval")
  write.csv(overall_raster,"Raster.csv")
  
  #Overall output compilation
  #Converting select variation prameters to RSD
  cvrsd <- function(data, sdindex, meanindex) {
    out = data
    count = 1
    for (i in 1:length(sdindex)) {
      out[,sdindex[count]] = data[,sdindex[count]]/data[,meanindex[count]]
      count = count + 1
    }
    return(out)
  }
  RSD <- c(4,6,12,13,14,28,32,34) #(BTB Var, DD Var, R30/60/90 Var, CTD Var, DS Var, Intensity Var)
  mean <- c(3,5,9,10,11,27,31,33) #(Columns of Means)
  overall_output <- cbind(overall_output,overall_output[,RSD])
  overall_output <- cvrsd(overall_output, RSD, mean)
  
  colnames(overall_output) <- c(Classifier, "Index",
                                "Interbeat_M", "Interbeat_V", "Depolarization_M", "Depolarization_V", "DepolarizationNCTD_M", "DepolarizationNCTD_V",
                                "Repolarization30_M", "Repolarization60_M", "Repolarization90_M", "Repolarization30_V", "Repolarization60_V", "Repolarization90_V", 
                                "Repolarization30NR60_M", "Repolarization30NR90_M", "Repolarization60NR90_M", "Repolarization30NR60_V", "Repolarization30NR90_V", "Repolarization60NR90_V", 
                                "Repolarization30NB_M", "Repolarization60NB_M", "Repolarization90NB_M", "Repolarization30NB_V", "Repolarization60NB_V", "Repolarization90NB_V", 
                                "CTD_M", "CTD_V", "CTDNB_M", "CTDNB_V", "Depolarization_Speed_M", "Depolarization_Speed_V",
                                "Intensity_M", "Intensity_V", "IntensityNCTD_M", "IntensityNCTD_V", "IntensityNB_M",
                                "CTD30_M", "CTD60_M", "CTD30_V", "CTD60_V", #Irrelevant parameters from this point on
                                "CTD30NCTD60_M", "CTD30NCTD90_M", "CTD60NCTD90_M", "CTD30NCTD60_V", "CTD30NCTD90_V", "CTD60NCTD90_V", 
                                "CTD30NB_M" , "CTD60NB_M", "CTD30NB_V", "CTD60NB_V",
                                "IntensityNN_M", "IntensityNBN_M", "IntensityNBN_V", "PeakIntensity_V", 
                                "Triangulation_M","Triangulation_V", "TriangulationNB_M", "TriangulationNB_V",
                                "Interbeat_SD","Depolarization_SD","R30_SD","R60_SD","R90_SD","CTD_SD","DS_SD","Intensity_SD" #SD of RSD converted parameters
                                )
  
  write.csv(overall_output,"Summary.csv")
  
  time = gsub(":","",substr(Sys.time(), 6,16))
  arrange = marrangeGrob(ROI_graph, nrow = 2, ncol = 2, top = textGrob(" ", gp = gpar(fontsize = 10)))
  ggsave(paste("Graphs",time,".pdf"), arrange, width = 24, height = 14, limitsize = F)
  dev.off()
  
}



################### Run this for Troubleshooting - Finding out which ROI has problems#########
# # # Troubleshooting
plot(dfc)
names(df[roi]) #Whatever the output is, minus 2, that's the ROI with an issue
defined_reads[sample]
head(dfc,4)
################### Run this for Troubleshooting - FInding out if peaks or first order derivative causes problems ###################
de_debug_1 #original
de_debug_2 #after doubled minima removal
de_outs #after no peak removal
outs_debug_1 #original
outs_debug_2 #after gstrict
outs_debug_3 #after no minima removal
outs #after no 1:1 minima removal


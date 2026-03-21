library(tidyverse)
library(gt)
library(openxlsx)
library(patchwork)


# read in sig treat result df
res <- read.xlsx('/Users/emilnetz/Desktop/cnv_sigs/MH/MH_therapy_outcome_pred/MH_allres_SigMedEnt_wCoGroups_20260220_2_polished.xlsx')

# select first 7 rows of res df 
res2 <- res %>% filter( q < 0.15)


# rename res2 columns 
res_plot <- res2 %>%
    rename(
        HR_log = coef,
        HR = `exp(coef)`,
        lower_CI_log = `coef.lower.95%`,
        upper_CI_log = `coef.upper.95%`,
        lower_CI = `exp(coef).lower.95%`,
        upper_CI = `exp(coef).upper.95%`,
    )


# plot main body of forest plot 
res_plot2 <- res_plot %>%
    mutate(across(
        c(HR, lower_CI, upper_CI), 
        ~ sprintf("%.2f", .x))) %>%
    mutate(HR_lab = paste0(HR, '(', lower_CI, '-', upper_CI,')')) %>% 
    mutate(q = case_when(
        q < 0.001 ~ '<0.001', 
        round(q, 2) == 0.05 ~ as.character(round(q , 3)), 
        q  < 0.01 ~ str_pad(
            as.character(round(q ,3)), 
            width = 4,
            pad = '0', 
            side = 'right'
        ),
        TRUE ~ str_pad(
            as.character(round(q ,2)),
            width = 4, 
            pad= '0', 
            side = 'right'
        )
    )) %>%
    mutate(total_cases = paste(total_cases)) %>%
    bind_rows(
        data.frame(
            SigMedEnt = "SigMedEnt",
            Signature = "Signature",
            Therapy = "Treatment", 
            Entity = "Entity",
            total_cases = "No. of Patients",
            HR_lab = "Hazard Ratio (95% CI)",
            lower_CI = "",
            upper_CI = "",
            q = "q-value"
        ) 
    ) %>%
    mutate(SigMedEnt = factor(SigMedEnt, levels = SigMedEnt)) %>% 
    mutate(SigMedEnt = fct_relevel(SigMedEnt, 'SigMedEnt')) %>%
    mutate(row_id = as.numeric(SigMedEnt))



# create forest plot 
p <- res_plot2 %>%
    ggplot(aes(y = fct_rev(SigMedEnt))) +
    theme_classic()+
    geom_rect(
    data = res_plot2 %>% filter(row_id %in% c(10, 8, 6, 4, 2)),
    aes(ymin = row_id - 0.5, ymax = row_id + 0.5),
    xmin = -Inf, xmax = Inf,
    fill = "grey90", inherit.aes = FALSE
  ) +
    geom_point(aes(x = HR_log), shape = 15, size = 3) + 
    geom_linerange(aes(xmin = lower_CI_log, xmax = upper_CI_log))+ 
    geom_vline(xintercept = 0, linetype = 'dashed') + 
    labs(x = 'Log Hazard Ratio', y = 'Signature/Treatment Interaction') + 
    coord_cartesian (ylim = c(1, 11), xlim = c(-0.08, 0.08)) + 
    annotate('text', x = -.045, y = 11, label = 'beneficial')+
    annotate('text', x = .045, y = 11, label = 'harmful')+ 
    theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.y= element_blank())


# plot left side of plot: Signature/Treatment + HR with CI 
p_left <- res_plot2 %>%
    ggplot(aes(y = fct_rev(SigMedEnt))) +
    geom_rect(
    data = res_plot2 %>% filter(row_id %in% c(10, 8, 6, 4, 2)),
    aes(ymin = row_id - 0.5, ymax = row_id + 0.5),
    xmin = -Inf, xmax = Inf,
    fill = "grey90", inherit.aes = FALSE
    )+
    # add signature treatment combination
    geom_text(aes(x = 0, label = Signature), hjust = 0, fontface= ifelse(res_plot2$Signature == "Signature", 'bold', 'plain'))+
    geom_text(aes(x = 0.6, label = Therapy), hjust = 0, fontface = ifelse(res_plot2$Signature == "Signature", 'bold', 'plain'))+
    geom_text(aes(x = 2.2, label = Entity), hjust = 0, fontface = ifelse(res_plot2$Signature == "Signature", 'bold', 'plain'))+
    geom_text(aes(x = 2.8, label= total_cases, hjust = 0, fontface = ifelse(res_plot2$Signature == "Signature", 'bold', 'plain')))+
    # add 
    geom_text(
        aes(x = 3.3, label = HR_lab), hjust = 0, fontface = ifelse(res_plot2$HR_lab == "Hazard Ratio (95% CI)", 'bold', 'plain'))+
    theme_void()+
    coord_cartesian(xlim = c(0,4))

# plot right part of plot: p values 
p_right <- res_plot2 %>% 
    ggplot()+
    geom_rect(
    data = res_plot2 %>% filter(row_id %in% c(10, 8, 6, 4, 2)),
    aes(ymin = row_id - 0.5, ymax = row_id + 0.5),
    xmin = -Inf, xmax = Inf,
    fill = "grey90", inherit.aes = FALSE
    )+
    geom_text(
    aes(x = 0, y = fct_rev(SigMedEnt), label = q),
    hjust = 0,
    fontface = ifelse(res_plot2$q == "q-value", "bold", "plain")
  ) +
  theme_void() 


# merge the three plots 

layout <- c(
  area(t = 0, l = 0, b = 30, r = 19), # r = 18
  area(t = 1, l = 18, b = 30, r = 27), # r = 9
  area(t = 0, l = 26, b = 30, r = 29) # r = 2
)

# plot rearrangement 

p_left + p + p_right + plot_layout(design = layout)


ggsave("/Users/emilnetz/Desktop/cnv_sigs/MH/figures/MH_treatpred_SigMedEnt_forestplot20260220.pdf", width=11.7, height=4)
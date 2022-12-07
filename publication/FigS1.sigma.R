library(tidyverse)

load("results/sleuth/sleuth_model_results/four_condition_contrasts/sleuth_mod_fourCondition.RData")

# Plot sigma-squared +/- co-var of interest
## Full model
fit_full <- sleuth_mod_fourCondition$fits$full_batch$summary%>%
  mutate(fit1 = "Full model") %>%
  select(target_id, sigma_sq, fit1) %>% 
  rename(sigma_sq1 = sigma_sq)

## Minus 1 covar
fit_all <- sleuth_mod_fourCondition$fits$full_noBatch$summary%>%
  mutate(fit2 = "Without batch") %>%
  select(target_id, sigma_sq, fit2) %>% 
  bind_rows(
    sleuth_mod_fourCondition$fits$full_noAge$summary%>%
      mutate(fit2 = "Without age") %>%
      select(target_id, sigma_sq, fit2)) %>% 
  bind_rows(
    sleuth_mod_fourCondition$fits$full_noSex$summary%>%
      mutate(fit2 = "Without sex") %>%
      select(target_id, sigma_sq, fit2)) %>% 
  rename(sigma_sq2=sigma_sq) %>% 
  full_join(fit_full) %>% 
  mutate(preferred_fit = case_when(sigma_sq1 < sigma_sq2 ~ fit1,
                                   sigma_sq2 < sigma_sq1 ~ fit2),
         preferred_fit = factor(preferred_fit, 
                                levels=c("Full model","Without age",
                                         "Without sex","Without batch")))

#### Summary table ####
fit_all %>% 
  group_by(fit2) %>% 
  count(preferred_fit) %>% 
  mutate(pct_without_covar = n[preferred_fit=="Full model"]/sum(n)*100)
  
## Add to facet labels
fit_all_lab <- fit_all %>% 
  mutate(facet_lab = case_when(fit2 == "Without age" ~ "28% of transcripts better fit\nwith age in model",
                               fit2 == "Without batch" ~ "81% of transcripts better fit\nwith batch in model",
                               fit2 == "Without sex" ~ "39% of transcripts better fit\nwith sex in model"))

#### Plot ####
fit_plot<-
  fit_all_lab%>%
  arrange(preferred_fit)%>%
  ggplot(aes(x=sigma_sq2, y = sigma_sq1))+
  geom_abline(slope = 1, intercept = 0, linetype="dashed")+
  geom_point(aes(color = preferred_fit))+
  labs(y = expression(paste(sigma^{2}, " full model" )),
       x = expression(paste(sigma^{2}, " model without co-variate")),
       color = "Preferred fit")+
  theme_bw()+
  coord_fixed() +
  scale_y_continuous(breaks = c(0,5,10))+
  facet_wrap(~facet_lab, nrow=1)
fit_plot

#### Save ####
ggsave("publication/FigS1.sigma.png", fit_plot, width=9, height=3)
ggsave("publication/FigS1.sigma.tiff", fit_plot, width=9, height=3)

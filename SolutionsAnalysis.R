# Project: Interlab ephys variance modeling
# Author: Dmitry Tebaykin
# Contact: dmitry.tebaykin@gmail.com
# Version: v1 (2 August 2017)

### Solution components analysis
### Requirements: DataProcessing.R

### Plot mean+/-SEM RMP for CA1 Pyr and Neostriatum medium spiny neuron
plot_data <- subset(data_models, NeuronName %in% c("Neostriatum medium spiny neuron","Hippocampus CA1 pyramidal cell"))
plot_data$NeuronName = as.factor(plot_data$NeuronName)
plot_data <- plot_data[!is.na(plot_data$rmp),]
plot_data <- plot_data[!is.na(plot_data$rmp_err),]
plot_data <- plot_data[with(plot_data, order(NeuronName, -rmp)),]
plot_data$index <- 1:nrow(plot_data)
plot_data$NeuronName = droplevels(plot_data$NeuronName)
levels(plot_data$NeuronName) = c("Hipp. CA1 pyr. cells", "Medium spiny neurons")
p = ggplot(plot_data, aes(x = index, y = rmp)) + 
  geom_point(aes(colour = NeuronName),  size = 0.3) +
  geom_pointrange(aes(ymin = rmp - rmp_err, ymax = rmp + rmp_err, colour = NeuronName), size = 0.1) +
  scale_color_manual(name = "Neuron type:", values=c("blue3", "darkgreen"),  guide='legend') +
  labs(x = "", y = "Vrest, mV") +
  theme_few(15) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
print(p)
ggsave(plot=p,height=4,width=6,dpi=200, filename= paste0(getwd(), "/Plots/exampleRMPnoMart.pdf"), useDingbats=FALSE)
rm(plot_data)

# Plots for external and internal solution concentrations
plot_data = subset(data_models, !duplicated(Pmid))
plot_data = plot_data[,relevant_solns]
plot_data = subset(plot_data, internal_0_EGTA < 50)
plot_data = subset(plot_data, internal_0_K > 100)
plot_data = subset(plot_data, external_0_Na > 100)
plot_data = subset(plot_data, external_0_Cl > 50)

summary(plot_data$external_0_Na)
summary(plot_data$external_0_Cl)
summary(plot_data$external_0_K)
summary(plot_data$external_0_Ca)
summary(plot_data$external_0_Mg)
summary(plot_data$external_0_glucose)

summary(plot_data$internal_0_Na)
summary(plot_data$internal_0_Cl)
summary(plot_data$internal_0_K)
summary(plot_data$internal_0_Ca)
summary(plot_data$internal_0_Mg)
summary(plot_data$internal_0_HEPES)
summary(plot_data$internal_0_EGTA)
summary(plot_data$internal_0_ATP)
summary(plot_data$internal_0_GTP)

colnames(plot_data) = new_relevant_solns
plot_data = melt(plot_data)
colnames(plot_data) = c("Compounds", "Concs")

# External and Internal compounds worth showing
plot_data = subset(plot_data, Concs < 200)

g = ggplot(subset(plot_data, Compounds %in% c("External [Na]","External [K]","External [Cl]","External [Ca]","External [Mg]","External [HEPES]", "External [glucose]") )) + 
  geom_histogram(aes(x = Concs, fill = Compounds), alpha = 0.5, binwidth = 0.5, position = "identity") +
  xlim(0, 15) +
  labs(x = "Compound concentration, mM", y = "Article count") +
  scale_fill_manual(values = c(brewer.pal(9, "Set1")[c(2:5,7)], "#00e6e6")) +
  theme_bw(20)
print(g)

my_pal = brewer.pal(9, "Set1")
my_pal[6] = my_pal[7]
my_pal[7] = "#00e6e6"
my_pal[9] = "#000000"
my_pal = c(my_pal, "#eef442", "#636363", "#cd93ff")
g = ggplot(subset(plot_data, Compounds %in% c("Internal [Na]","Internal [K]","Internal [Cl]","Internal [Ca]","Internal [Mg]", "Internal [HEPES]","Internal [EGTA]","Internal [ATP]","Internal [GTP]", "Internal [BAPTA]", "Internal [Cs]", "Internal [glucose]") )) + 
  geom_histogram(aes(x = Concs, fill = Compounds), alpha = 0.5, binwidth = 1, position = "identity") +
  labs(x = "Compound concentration, mM", y = "Article count") +
  scale_fill_manual(values = my_pal) +
  theme_bw(20)
print(g)

g = ggplot(subset(plot_data, Compounds %in% c("Internal [Na]","Internal [K]","Internal [Cl]","Internal [Ca]","Internal [Mg]", "Internal [HEPES]","Internal [EGTA]","Internal [ATP]","Internal [GTP]", "Internal [BAPTA]", "Internal [Cs]", "Internal [glucose]") )) + 
  geom_histogram(aes(x = Concs, fill = Compounds), alpha = 0.5, binwidth = 0.5, position = "identity") +
  xlim(0, 15) +
  labs(x = "Compound concentration, mM", y = "Article count") +
  scale_fill_manual(values = my_pal) +
  theme_bw(20)
print(g)

# Multiplot of solution concs
p=ggplot(plot_data) +
  geom_histogram(aes(x = Concs, fill = Compounds), alpha = 0.5, bins = 50, position = "identity") +
  labs(x = "Compound concentration, mM", y = "Article count") +
  theme_bw(15) +
  facet_wrap(~Compounds, scales = "free") +
  theme(legend.position="none")
print(p)
ggsave(plot=p,height=10,width=10,dpi=200, filename= paste0(getwd(), "/Plots/SolnConcs.pdf"), useDingbats=FALSE)

# Internal Na vs year published 
plot_data = subset(data_models, !duplicated(Pmid))
p = ggplot(plot_data, aes(x = PubYear, y = internal_0_Na, group = PubYear)) + 
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(alpha = 0.3, position = position_jitter(width = 0.2)) +
  
  geom_pointrange(stat = "summary", color = "darkred", alpha = 0.6,
                  fun.ymin = mean, #function(z) {quantile(z,0.25)},
                  fun.ymax = mean, #function(z) {quantile(z,0.75)},
                  fun.y = mean) +
  labs(y = "Na+int", x = "Publication year") +
  scale_x_continuous(breaks = seq(min(plot_data$PubYear), max(plot_data$PubYear), by = 1)) +
  theme_few(15)
print(p)
ggsave(plot=p,height=4,width=8,dpi=200, filename= paste0(getwd(), "/Plots/NaVsYear.pdf"), useDingbats=FALSE)

# Unique solutions exploration
table(duplicated(plot_data[,c("external_0_Na", "external_0_Cl", "external_0_K", "external_0_Mg", "external_0_Ca")]))
table(duplicated(plot_data[,c("internal_0_Na", "internal_0_Cl", "internal_0_K", "internal_0_Mg", "internal_0_Ca")]))

table(duplicated(plot_data[,c("internal_0_Na", "internal_0_Cl", "internal_0_K", "internal_0_Mg", "internal_0_Ca",
                              "external_0_Na", "external_0_Cl", "external_0_K", "external_0_Mg", "external_0_Ca")]))
library(ggplot2)
library(dplyr)

### Pie chart

df <- data.frame(dataset = 'Ureter10_Urothelial',
                 category = c('Exonic','Flank','Intronic','TSS-proximal','UTR3'),
                 count = c(6029,5468,2087,1651,16992))
cols <- c('UTR3' = '#762a83','TSS-proximal'='#74add1','Exonic'='#ffffbf',
          'Intronic'='#f46d43', 'Flank' = '#a50026')

df <- df %>% group_by(dataset) %>%
  mutate(percentage = (count / sum(count)) * 100)

df$percentage <- paste0(round(df$percentage, digits = 1),'%')

df$category <- factor(df$category, levels = c('Flank','Intronic','Exonic','TSS-proximal','UTR3'))

# pie chart
p <- ggplot(df, aes(x="", y=count, fill=category)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values=cols) +
  theme_void() + geom_text(aes(label = percentage),
                           position = position_stack(vjust = 0.5))

p <- ggplot(df, aes(x="", y=count, fill=category)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values=cols) +
  theme_void()

p <- p + theme(text = element_text(size = 14)) 

ggsave('~/Desktop/Ureter10_scPASU_run/outputs/urothelial/3h_fragmented_peaks_to_merge/Ureter10_urothelial_final_peaks_piechart_nolabel.png',
       width = 8, height = 8, units = 'in', bg = 'white', p)

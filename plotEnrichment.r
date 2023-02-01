plotEnrichment <-function (data = data, type = c("CGI", "Genic")) {
cpgPal <- c("#D43F3AFF","#EEA236FF","#5CB85CFF","#46B8DAFF")
genicPal <- c("#800000FF","#767676FF","#FFA319FF","#8A9045FF","#155F83FF","#C16622FF","#250F20FF")
stopifnot(type %in% c("CGI", "Genic"))
data <- data %>% mutate(OR = case_when(OR < 1 ~ -1/OR, OR >= 1 ~ OR)) %>% dplyr::mutate(signif = dplyr::case_when(fdr <= 0.05 ~ 1, fdr > 0.05 ~ 0))
p <- ggplot(data = data, aes(x = Annotation, y = OR, fill = Annotation)) + geom_bar(stat = "identity", color = "Black") + coord_flip() + labs(y = "Fold Enrichment", x = element_blank()) + theme_classic() + theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20), strip.text = element_text(size = 20), legend.text = element_text(size = 20), legend.position = "none") + scale_y_continuous(expand = c(0.1, 0.1)) + scale_x_discrete(limits = data$Annotation %>% as.factor() %>% levels() %>% rev()) + geom_hline(yintercept = 0) + geom_text(data = data[(data$signif == 1 & data$OR > 0),], label = "*", size = 12, show.legend = FALSE, nudge_y = 0.5, nudge_x = -0.09) + geom_text(data = data[(data$signif == 1 & data$OR < 0), ], label = "*", size = 8, show.legend = FALSE, nudge_y = -0.5, nudge_x = -0.09)
if (type == "CGI") {
p <- p + scale_fill_manual(values = c(cpgPal), breaks = data$Annotation %>% as_factor() %>% levels(), name = "Annotation")
}
else if (type == "Genic") {
p <- p + scale_fill_manual(values = genicPal)
}
return(p)
}

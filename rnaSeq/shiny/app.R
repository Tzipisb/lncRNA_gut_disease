library(shiny)
library(shinythemes)
library(ggplot2)
library(shinysky) 
library(patchwork)
source('helpers.R')


data_file = 'data/main6_lncRNA_TPM_v8.rds'
tpm = readRDS(data_file)
names(tpm)[names(tpm) == 'Xist'] = 'Xist.exon'
# names(tpm) = lapply(names(tpm), tolower)
# names(tpm) = lapply(names(tpm), make.names)

gene_names <- gsub("\\.", "-", names(tpm)[8:length(names(tpm))]) 

# Define UI ----
ui =  fluidPage(theme = shinytheme("flatly"),
  titlePanel("LncRNA Atlas of Inflammatory Diseases Along the GI Tract"),
  sidebarLayout(
    sidebarPanel(h3("lncRNA name"),  
                 # textInput("lncRNA", label = h3("lncRNA name"), value = "GATA6-AS1"),
                 textInput.typeahead(  id = "lncRNA"
                                       , valueKey = "lncRNA"
                                       , local = data.frame(lncRNA = c(gene_names))
                                       , tokens = c(1:length(gene_names))
                                       , placeholder = "GATA6-AS1"
                                       , template = HTML("<p>{{lncRNA}}</p>")),
                 
                 hr(),
                 # fluidRow(column(3, verbatimTextOutput("lncRNA"))),
                 
                 radioButtons("COHORT", 
                              # label = h3("Disease and gut location"),
                              label = h3("Data Type"),
                              choices = list("Ulcerative Colitis (UC), Rectum" = "PROTECT", 
                                             "Crohn's Disease (CD), Ileum" = "SOURCE", 
                                             "Celiac, Duodenum" = "SEEM", 
                                             "GI tract expresion\nin controls" = 'Controls'), 
                              selected = "PROTECT"),
                 
                 hr(),
                 # fluidRow(column(3, verbatimTextOutput("Cohort")))
                 ),
    mainPanel("", 
              h3(""),  
              # textOutput("cohort_text"), 
              br(),
              # textOutput("lncRNA_text"), 
              plotOutput("boxplot", width = "100%"),
              # plotOutput("boxplot2", width = "100%"),
              # fluidRow(splitLayout(cellWidths = c("50%", "50%"), plotOutput("boxplot"), plotOutput("boxplot2")) ), 
              textOutput("stats_text"), 
              htmlOutput("cohort_full_text"), 
              ),
  )
)

# Define server logic ----
server = function(input, output) 
{
  output$cohort_text <- renderText({ sprintf("Cohort: %s",input$COHORT) })
  output$cohort_full_text <- renderUI({ if (input$COHORT == 'PROTECT')
    {
    HTML(paste('',
          'PROTECT is a multicenter inception cohort study including treatment naïve children aged 4–17 years with a diagnosis of UC. Rectal mucosal biopsies from a representative sub-cohort of 206 PROTECT UC patients and 20 age and gender matched non-IBD controls, underwent high coverage transcriptomic profiling using Illumina RNAseq.', 
          '',
          'RISK cohort is a treatment naïve pediatric IBD cohort. 55 control and 43 UC rectal biopsies of this cohort were used here.',
          '',
          'Reference:',
          'Hyams, J. S. et al. Clinical and biological predictors of response to standardised paediatric colitis therapy (PROTECT): a multicentre inception cohort study. Lancet 393, 1708-1720, doi:10.1016/S0140-6736(18)32592-3 (2019).',
          'Haberman, Y. et al. Ulcerative colitis mucosal transcriptomes reveal mitochondriopathy and personalized mechanisms underlying disease severity and treatment response. Nat Commun 10, 38, doi:10.1038/s41467-018-07841-3 (2019).',
          sep="<br/>"))
    } else if ( input$COHORT == 'SEEM' )
    {
      HTML(paste('',
                 'SEEM celiac cohort included 17 celiac subjects and 25 controls from Cincinnati Children’s Hospital Medical Center (CCHMC). Controls were subjects who were investigated for various gastrointestinal symptoms including abdominal pain, but had normal endoscopic and histologic findings. Celiac disease diagnosis was based on previously described algorithms including positive IgA autoantibodies against tissue transglutaminase (anti-TTG) and histologic features.',
                 '',
                 'PRJNA52875 is a published celiac cohort. 15 control and 12 activae celiac duodenum biopsies of this cohort were used here.',
                 '',
                 'Reference:',
                 'Haberman, Y. et al. Mucosal Genomics Implicate Lymphocyte Activation and Lipid Metabolism in Refractory Environmental Enteric Dysfunction. Gastroenterology 160, 2055-2071 e2050, doi:10.1053/j.gastro.2021.01.221 (2021).',
                 'Leonard, M. M. et al. RNA sequencing of intestinal mucosa reveals novel pathways functionally linked to celiac disease pathogenesis. PLoS One 14, e0215132, doi:10.1371/journal.pone.0215132 (2019).',
                  sep="<br/>"))
    } else if ( input$COHORT == 'SOURCE' )
    {
      HTML(paste('',
                'SOURCE cohort included 18 treatment naïve Crohn disease and 25 non-IBD controls from the Sheba Medical Center.',
                '',
                'RISK cohort is a treatment naïve pediatric IBD cohort. 47 control and 213 iCD ileal biopsies of this cohort were used here.',
                '',
                'Reference:',
                'Haberman, Y. et al. Pediatric Crohn disease patients exhibit specific ileal transcriptome and microbiome signature. J Clin Invest 124, 3617-3633, doi:10.1172/JCI75436 (2014).',
                sep="<br/>"))
      
    } else if ( input$COHORT == 'Controls' ) 
    {
      HTML(paste('',
                 'PROTECT is a multicenter inception cohort study including treatment naïve children aged 4–17 years with a diagnosis of UC. Rectal mucosal biopsies from a representative sub-cohort of 206 PROTECT UC patients and 20 age and gender matched non-IBD controls, underwent high coverage transcriptomic profiling using Illumina RNAseq.',
                 '',
                 'RISK cohort is a treatment naïve pediatric IBD cohort. 55 control and 43 UC rectal biopsies of this cohort were used here.',
                 '',
                 'Reference:',
                 'Hyams, J. S. et al. Clinical and biological predictors of response to standardised paediatric colitis therapy (PROTECT): a multicentre inception cohort study. Lancet 393, 1708-1720, doi:10.1016/S0140-6736(18)32592-3 (2019).',
                 'Haberman, Y. et al. Ulcerative colitis mucosal transcriptomes reveal mitochondriopathy and personalized mechanisms underlying disease severity and treatment response. Nat Commun 10, 38, doi:10.1038/s41467-018-07841-3 (2019).',
                 '',
                 'SOURCE cohort included 18 treatment naïve Crohn disease and 25 non-IBD controls from the Sheba Medical Center.',
                 '',
                 'RISK cohort is a treatment naïve pediatric IBD cohort. 47 control and 213 iCD ileal biopsies of this cohort were used here.',
                 '',
                 'Reference:',
                 'Haberman, Y. et al. Pediatric Crohn disease patients exhibit specific ileal transcriptome and microbiome signature. J Clin Invest 124, 3617-3633, doi:10.1172/JCI75436 (2014).',
                 '',
                 'SEEM celiac cohort included 17 celiac subjects and 25 controls from Cincinnati Children’s Hospital Medical Center (CCHMC). Controls were subjects who were investigated for various gastrointestinal symptoms including abdominal pain, but had normal endoscopic and histologic findings. Celiac disease diagnosis was based on previously described algorithms including positive IgA autoantibodies against tissue transglutaminase (anti-TTG) and histologic features.',
                 '',
                 'PRJNA52875 is a published celiac cohort. 15 control and 12 activae celiac duodenum biopsies of this cohort were used here.',
                 '',
                 'Reference:',
                 'Haberman, Y. et al. Mucosal Genomics Implicate Lymphocyte Activation and Lipid Metabolism in Refractory Environmental Enteric Dysfunction. Gastroenterology 160, 2055-2071 e2050, doi:10.1053/j.gastro.2021.01.221 (2021).',
                 'Leonard, M. M. et al. RNA sequencing of intestinal mucosa reveals novel pathways functionally linked to celiac disease pathogenesis. PLoS One 14, e0215132, doi:10.1371/journal.pone.0215132 (2019).',
                 sep="<br/>"))
      }}
      )
  output$lncRNA_text <- renderText({ sprintf("lncRNA: %s",input$lncRNA) })
  
  output$stats_text = renderText({ifelse(input$COHORT!='Controls',
                                         '*p<0.05, **p<0.01 ***p<0.001, calculated using Mann-Whitney test between cases and controls.',
                                         '')})
  
  tpm2 = reactive( { tpm[tpm$Cohort == input$COHORT,] } )
  tpm_control = reactive( { tpm[tpm$Dx == 'Control',] } )
  lnc_pos = reactive({ 
    lnc = make.names(input$lncRNA)
    grep(lnc,names(tpm),ignore.case = T ) 
    })

  tpm2_val = reactive( { temp = input$COHORT
                         temp = gsub('PROTECT','RISK Rectal',temp)
                         temp = gsub('SOURCE','RISK ileal',temp)
                         temp = gsub('SEEM','PRJNA52875',temp)
                         tpm[tpm$Cohort == temp,] } )
  lnc_pos = reactive({ 
    lnc = make.names(input$lncRNA)
    grep(lnc,names(tpm),ignore.case = T ) 
  })
  # output$cohort_text <- renderText({ sprintf("Cohort: %s",names(tpm[300])) })
  # output$lncRNA_text <- renderText({ sprintf("lncRNA: %s",lnc_pos() ) })
  
  output$boxplot = renderPlot( { # ggplot()
    lnc = lnc_pos()
    if(input$lncRNA == '')
    {
      lnc=grep( make.names('GATA6-AS1'),names(tpm),ignore.case = T ) 
    }
    if (input$COHORT != 'Controls')
    {
      
      g1 = make_disease_boxplot(tpm2(), lnc)
      g2 = make_disease_boxplot(tpm2_val(), lnc)
      return(g1 + (g2+  ylab('')) )
    } else
    {
      tpmc = tpm_control()
      tpmc$Source = factor(tpmc$Source, levels = c('Rectum','Ileum','Duodenum'))
      tpmc$Cohort = factor(tpmc$Cohort, levels = c('PROTECT','SOURCE','SEEM','RISK Rectal','RISK ileal','PRJNA52875'))
      g = ggplot(tpmc, aes(x=Cohort, y=tpmc[,lnc])) + geom_boxplot(fill = 'gray70') + 
        ylab(sprintf('%s (TPM) in\nControls',names(tpmc)[lnc] )) + 
        theme_bw() + facet_wrap(~Source, scales = 'free_x') + 
        theme(axis.title = element_text(size=18) , 
              axis.text.x = element_text(size=14, angle = 45, hjust = 1), 
              axis.text.y = element_text(size=14), 
              strip.text.x = element_text(size = 14))
      return(g)
    }
  }, height = 250, width = 500 )
 
  
}

# Run the app ----
shinyApp(ui = ui, server = server)

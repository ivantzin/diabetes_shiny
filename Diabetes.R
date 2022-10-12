list.of.packages <- c("deSolve","tidyr", "shiny", "dplyr", "ggplot2", "cowplot", "shinydashboard","ggrepel","gridExtra","reshape2","readxl","forcats","scales")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
invisible(lapply(list.of.packages, library, character.only = TRUE))

#import parameters
coeffs <- read_excel("Parameters/parameters_feb.xlsx")
coeffs$Description<-NULL
coeffs$Type<-NULL
coeffs <- coeffs %>% spread(Rate,Value)

#Model equations----
diab <- function(time, state, parameters)  {
  with(as.list(c(state, parameters)), {
    
    #ty names timesteps initing from 2000
    ty=time+inityear
    
    #Population equations----
    
    #initial states At risk w/w.o No Risk Factor Reduction & prev. intervention
    dA= i3*S_PD_A*PD+i3*S_ND_A*ND+R_UND_A*S_UND_A*UND-(i1*S_A_SUS+i1*S_A_ND+(1-i1)*S_A_UND+(1-i1)*S_A_UD)*A
    dANRFR=(1-i3)*S_ND_ANRFR*ND+(1-i3)*S_PD_ANRFR*PD+ R_UND_ANRFR*S_UND_ANRFR*UND -(i1*S_ANRFR_SUS+i1*S_ANRFR_ND+(1-i1)*S_ANRFR_UND+(1-i1)*S_ANRFR_UD)*ANRFR
    #first screening
    dSUS=i1*S_A_SUS*A+i1*S_ANRFR_SUS*ANRFR-(S_ND+S_PD+S_D1+S_D2)*i2*SUS
    #second screening
    dND=i1*S_A_ND*A+i1*S_ANRFR_ND*ANRFR+i2*S_ND*SUS-(i3*S_ND_A+(1-i3)*S_ND_ANRFR)*ND
    dPD=i2*S_PD*SUS-(i3*S_PD_A+(1-i3)*S_PD_ANRFR)*PD
    dD1=i2*S_D1*SUS-(R_D1HC*S_D1HC+R_D1LC*S_D1LC+R_D1U*S_D1U)*D1
    dD2=i2*S_D2*SUS-(R_D2HC*S_D2HC+R_D2LC*S_D2LC+R_D2U*S_D2U)*D2
    
    #undiagnosed population
    dUND=(1-i1)*S_A_UND*A+(1-i1)*S_ANRFR_UND*ANRFR-(R_UND_A*S_UND_A+R_UND_ANRFR*S_UND_ANRFR + R_UND_UD*S_UND_UD )*UND
    dUD=(1-i1)*S_A_UD*A+(1-i1)*S_ANRFR_UD*ANRFR + R_UND_UD*S_UND_UD*UND - R_UD_M*S_UD_M*UD - (R_UD_DR*S_UD_DR+R_UD_DNEUR*S_UD_DNEUR+R_UD_DNEPH*S_UD_DNEPH+R_UD_CVD*S_UD_CVD)*UD
    
    #D1 states 
    dD1HC=R_D1HC*S_D1HC*D1
    dD1LC=R_D1LC*S_D1LC*D1-R_D1LC_D2LC*S_D1LC_D2LC*D1LC
    dD1U=R_D1U*S_D1U*D1-R_D1U_D2U*S_D1U_D2U*D1U
    
    #D2 states
    dD2HC=R_D2HC*S_D2HC*D2
    dD2LC=R_D2LC*S_D2LC*D2+R_D1LC_D2LC*S_D1LC_D2LC*D1LC-R_D2LC_M*S_D2LC_M*D2LC
    dD2U=R_D1U_D2U*S_D1U_D2U*D1U+R_D2U*S_D2U*D2-(R_D2U_M*S_D2U_M+R_D2U_CVD*S_D2U_CVD+R_D2U_DNEPH*S_D2U_DNEPH+R_D2U_DNEUR*S_D2U_DNEUR+R_D2U_DR*S_D2U_DR)*D2U
    
    
    #Monitored state
    
    dM=R_UD_M*S_UD_M*UD+R_D2U_M*S_D2U_M*D2U+R_D2LC_M*S_D2LC_M*D2LC-(R_M_CVD*S_M_CVD+R_M_DNEPH*S_M_DNEPH+R_M_DNEUR*S_M_DNEUR+R_M_DR*S_M_DR)*M
    
    #Complications
    
    dM_CVD=R_M_CVD*S_M_CVD*M  #-mu1_CVD*UD_CVD-mu*M_CVD
    dM_DNEPH=R_M_DNEPH*S_M_DNEPH*M #-mu1_DNEPH*UD_DNEPH-mu*M_DNEPH
    dM_DNEUR=R_M_DNEUR*S_M_DNEUR*M  #-mu1_DNEUR*UD_DNEUR-mu*M_DNEUR
    dM_DR=R_M_DR*S_M_DR*M  #-mu1_DR*UD_DR-mu*M_DR
    
    dUD_CVD= R_D2U_CVD*S_D2U_CVD*D2U + R_UD_CVD*S_UD_CVD*UD #- mu2_CVD*UD_CVD - mu*M_CVD
    dUD_DNEPH=R_D2U_DNEPH*S_D2U_DNEPH*D2U + R_UD_DNEPH*S_UD_DNEPH*UD #- mu2_DNEPH*UD_DNEPH - mu*M_DNEPH
    dUD_DNEUR= R_D2U_DNEUR*S_D2U_DNEUR*D2U + R_UD_DNEUR*S_UD_DNEUR*UD #- mu2_DNEUR*UD_DNEUR - mu*M_DNEUR
    dUD_DR= R_D2U_DR*S_D2U_DR*D2U + R_UD_DR*S_UD_DR*UD #- mu2_DR*UD_DR - mu*M_DR
    
    
    #Model output----
    return(list(c(dA,dANRFR,dSUS,dND, dPD,dD1, dD2, dUND, dUD,dD1HC,dD1LC,dD1U,dD2HC,dD2LC,dD2U,dM,dM_CVD,dM_DNEPH,dM_DNEUR,dM_DR,dUD_CVD,dUD_DNEPH,dUD_DNEUR,dUD_DR)))
  })
}


ui <- dashboardPage(
  dashboardHeader(disable = TRUE),
  dashboardSidebar(
    sliderInput("i1",
                "Screening coverage (%) ",
                min = 0, max = 100, value = 50),
    sliderInput("i3",
                "Risk factor reduction intervention coverage (%)",
                min = 0, max = 100, value = 20),
    sliderInput("treat",
                "Treated patients with Diabetes (%)",
                min = 0, max = 100, value = 50),
 #   sliderInput("treat",
  #              "Treated patients with D2 (%)",
   #             min = 0, max = 100, value = 20),
    sliderInput("hcD",
                "Patients under High cost treatment  (%)",
                min = 0, max = 100, value = 20),
    sliderInput("monit",
                "Monitored patients (%)",
                min = 0, max = 100, value = 20)
  ),
  dashboardBody(
    tags$head(tags$style(HTML('
                              /* body */
                              .content-wrapper, .right-side {
                              background-color: #fffff8;
                              }                              
                              '))),
    
    #    mainPanel(
    fluidRow(plotOutput("healthplot")),
    fluidRow(plotOutput("econplot")),
    
    br(),
    verbatimTextOutput(""),
    fluidRow(
      # Dynamic valueBoxes
      valueBoxOutput("Profit", width = 4)
    )
  )
)
#
# Define server 
#
server <- function(input, output) {
  # Create reactive input
  dataInput <- reactive({
    init       <-  c(
      #initial conditions
      A=coeffs$A,	ANRFR=coeffs$ANRFR,	SUS=coeffs$SUS,	ND=coeffs$ND,	PD=coeffs$PD,	D1=coeffs$D1,	D2=coeffs$D2,	UND=coeffs$UND,	UD=coeffs$UD,	D1HC=coeffs$D1HC,	D1LC=coeffs$D1LC,	D1U=coeffs$D1U,	D2HC=coeffs$D2HC,	D2LC=coeffs$D2LC,	D2U=coeffs$D2U,	M=coeffs$M,	M_CVD=coeffs$M_CVD,	M_DNEPH=coeffs$M_DNEPH,	M_DNEUR=coeffs$M_DNEUR,	M_DR=coeffs$M_DR,	UD_CVD=coeffs$UD_CVD,	UD_DNEPH=coeffs$UD_DNEPH,	UD_DNEUR=coeffs$UD_DNEUR,	UD_DR=coeffs$UD_DR)
    
    parameters <-c(
      #interventions 
      i1=input$i1/100, i3=input$i3/100, i2=coeffs$i2,	
      #rates
      R_D1HC=coeffs$R_D1HC,	R_D1LC=coeffs$R_D1LC,	R_D1LC_D2LC=coeffs$R_D1LC_D2LC,	R_D1U=coeffs$R_D1U,	R_D1U_D2U=coeffs$R_D1U_D2U,	R_D2HC=coeffs$R_D2HC,	R_D2LC=coeffs$R_D2LC,	R_D2U=coeffs$R_D2U,	R_D2LC_M=coeffs$R_D2LC_M,	R_D2U_M=coeffs$R_D2U_M,	R_UND_A=coeffs$R_UND_A,	R_UND_ANRFR=coeffs$R_UND_ANRFR,	R_UND_UD=coeffs$R_UND_UD,	R_M_CVD=coeffs$R_M_CVD,	R_M_DNEPH=coeffs$R_M_DNEPH,	R_M_DNEUR=coeffs$R_M_DNEUR,	R_M_DR=coeffs$R_M_DR,	R_D2U_CVD=coeffs$R_D2U_CVD,	R_D2U_DNEPH=coeffs$R_D2U_DNEPH,	R_D2U_DNEUR=coeffs$R_D2U_DNEUR,	R_D2U_DR=coeffs$R_D2U_DR,	R_UD_M=coeffs$R_UD_M,	R_UD_CVD=coeffs$R_UD_CVD,	R_UD_DNEPH=coeffs$R_UD_DNEPH,	R_UD_DNEUR=coeffs$R_UD_DNEUR,	R_UD_DR=coeffs$R_UD_DR,	
       #SHARES
                   S_A_SUS=coeffs$S_A_SUS,	S_A_ND=coeffs$S_A_ND,	S_A_UND=coeffs$S_A_UND,	S_A_UD=coeffs$S_A_UD,	S_ND=coeffs$S_ND,	S_PD=coeffs$S_PD,	S_D1=coeffs$S_D1,	S_D2=coeffs$S_D2,	S_ANRFR_SUS=coeffs$S_ANRFR_SUS,	S_ANRFR_ND=coeffs$S_ANRFR_ND,	S_ANRFR_UND=coeffs$S_ANRFR_UND,	S_ANRFR_UD=coeffs$S_ANRFR_UD,	S_D1LC_D2LC=coeffs$S_D1LC_D2LC,	S_D1U_D2U=coeffs$S_D1U_D2U,	S_D2U=coeffs$S_D2U,	S_D2HC=coeffs$S_D2HC,	S_D2LC=coeffs$S_D2LC,	S_ND_A=coeffs$S_ND_A,	S_PD_A=coeffs$S_PD_A,	S_ND_ANRFR=coeffs$S_ND_ANRFR,	S_PD_ANRFR=coeffs$S_PD_ANRFR,	S_M_CVD=coeffs$S_M_CVD,	S_M_DNEPH=coeffs$S_M_DNEPH,	S_M_DNEUR=coeffs$S_M_DNEUR,	S_M_DR=coeffs$S_M_DR,	S_UND_A=coeffs$S_UND_A,	S_UND_ANRFR=coeffs$S_UND_ANRFR,	S_UND_UD=coeffs$S_UND_UD,
      #shares based on user's slider inputs. 
      S_D1HC = input$hcD*input$treat/10000,	S_D1LC = (100-input$hcD)*input$treat/10000,	S_D1U = 1-(input$treat)/100,S_D2U =  1-(input$treat)/100,	S_D2HC =input$hcD*input$treat/10000,	S_D2LC = (100-input$hcD)*input$treat/10000,	S_D2LC_M =  input$monit/100, S_D2U_M = input$monit/100,	S_D2U_CVD = (1-input$monit/100)/4 ,	S_D2U_DNEPH =(1-input$monit/100)/4 ,	S_D2U_DNEUR = (1-input$monit/100)/4 ,	S_D2U_DR = (1-input$monit/100)/4,S_UD_M = input$monit/100,	S_UD_CVD = (1-input$monit/100)/4,	S_UD_DNEPH = (1-input$monit/100)/4,	S_UD_DNEUR = (1-input$monit/100)/4,	S_UD_DR = (1-input$monit/100)/4,inityear=2000)
    ## Time frame
    years<-30
    times <- seq(0, years, 1/365.25)
    #times<-seq(0, 40, 1/365.25)
    ## Solve using ode (General Solver for Ordinary Differential Equations)
    run_d<-ode(times=times, y=init, func=diab, parms=parameters)
    #    out
    out<-as.data.frame(run_d)
    out$Diabetes1<-out$D1+out$D1HC+out$D1LC+out$D1U
    out$Diabetes2<-out$D2+out$D2HC+out$D2LC+out$D2U
    out$NonDiab<-out$UND+out$ND+out$A+out$ANRFR
    out$PreDiab<-out$PD
    out$Undiagnosed_D<-out$UD
    out$Monitored<-out$M
    out$CompM<-out$M_CVD+out$M_DNEPH+out$M_DNEUR+out$M_DR
    out$CompUn<-out$UD_CVD+out$UD_DNEPH+out$UD_DNEUR+out$UD_DR
    out <- out[-c(2:25)]
    out[,-1] <-round(out[,-1],0) 
    out<-as.data.frame(out)
    out <- out %>% rowwise() %>%
      mutate(other = 10000000-sum(c_across(Diabetes1:CompUn)))
  })
  
  output$healthplot <- renderPlot({

    #Creating data for each dataset
    out2 <-
      dataInput()  %>% select(time,Diabetes1:Undiagnosed_D) %>%
      gather(State, Population, -time) %>%
      mutate(
        id = row_number(),
        key = recode(
          State,
          Diabetes1 = "Diabetes 1",
          Diabetes2 = "Diabetes 2",
          NonDiab= "Non-diabetic",
          PreDiab = "Prediabetic",
          Undiagnosed_D="Undiagnosed diabetic"
        )
      )
    #ggplot(out,aes(x=time,y=Population,colour=key,group=key)) + geom_line()
    out3 <-
      dataInput()  %>% select(time,Monitored:CompUn) %>%
      gather(State, Population, -time) %>%
      mutate(
        id = row_number(),
        key = recode(
          State,
          Monitored="Monitored D2 w/o comp.",
          CompM="Monitored D2 w/comp.", 
          CompUn="Unmonitored D2 w/comp." 
        )
      )

    ###Create graphs  
    
    p1<-ggplot(data = out2,
               aes(
                 x = time,
                 y = Population,
                 group = key,
                 col = key,
                 label = key,
                 data_id = id
               ))  + ylab("Population") + xlab("Time (years)") +
      geom_line(size = 1) +
      geom_text_repel(
        data = subset(out2, time == max(time)),
        size = 4,
        nudge_x = 0,
        segment.color = NA
      ) +
      theme(legend.position = "none") +
      #scale_colour_manual(values = c("Magenta", "orange", "blue", "violet", "green","black","purple","red")) +
      theme(
        rect=element_rect(size=0),
        legend.position="none",
        legend.key = element_rect(fill = "transparent", colour = "transparent")
      ) +ggtitle("Diabetic status of population")
    
    p2<- ggplot(data = out3,
                aes(
                  x = time,
                  y = Population,
                  group = key,
                  col = key,
                  label = key,
                  data_id = id
                ))  + ylab("Population") + xlab("Time (years)") +
      geom_line(size = 1) +
      geom_text_repel(
        data = subset(out3, time == max(time)),
        size = 4,
        nudge_x = 0,
        segment.color = NA
      ) +
      theme(legend.position = "none") +
      #scale_colour_manual(values = c("Magenta", "orange", "blue", "violet", "green","black","purple","red")) +
      theme(
        rect=element_rect(size=0),
        legend.position="none",
        legend.key = element_rect(fill = "transparent", colour = "transparent")
      ) +ggtitle("Associated health outcomes")
  
  
    
    #Arrangement of graphs
    grid.arrange(p1,p2,ncol=1)
  }) 
  
#second plot
  
  output$econplot <- renderPlot({
    
    #Creating data for each dataset

    piedata <-
      dataInput()  %>% select(Diabetes1:other) %>%   
      gather(population, value, c(Diabetes1:other))   %>%
      group_by(population) %>% do(tail(., n=1))
    piedata$value <- piedata$value/10000000
    #  piedata$value <-piedata %>%  mutate_each(funs(try(percent(.), silent=TRUE)), -population)
    
    temp <- dataInput() 
    n_last<-length(temp$other)
    
    econdata<-data.frame(
      #costs
        (input$i1/100)*coeffs$p_screen*temp$NonDiab,
        (input$i3/100)*coeffs$p_rfr*temp$NonDiab, 
coeffs$p_hct*(temp$Diabetes1+temp$Diabetes2)*(input$treat/100)*(input$hcD/100)+ coeffs$p_lct*(temp$Diabetes1+temp$Diabetes2)*(input$treat/100)*((1-input$hcD)/100)+ coeffs$p_hct*(temp$Monitored)*(input$treat/100)*(input$hcD/100)+ coeffs$p_lct*(temp$Monitored)*(input$treat/100)*((1-input$hcD)/100),
        coeffs$p_monit*temp$Monitored+coeffs$p_monit*temp$CompM,
      #benefits
      coeffs$p_yll*temp$NonDiab+coeffs$p_yll_pd*temp$PreDiab,
      coeffs$p_abs*temp$NonDiab,
      coeffs$p_prs*temp$NonDiab
        )
    colnames(econdata) <- c("Screening cost", "RFR cost", "Diabetes treatment cost","Monitoring cost","YLL savings", "Absenteeism cost savings", "Presenteeism cost savings")
    econdata$n <- 1:nrow(econdata)
    econdata$irate <- (1/(1+.05))^econdata$n
    econdata$`Screening cost`<- econdata$`Screening cost`/ econdata$irate
    econdata$`RFR cost`<- econdata$`RFR cost`/ econdata$irate
    econdata$`Diabetes treatment cost`<- econdata$`Diabetes treatment cost`/ econdata$irate
    econdata$`Monitoring cost`<- econdata$`Monitoring cost`/ econdata$irate
    econdata$`YLL savings`<- econdata$`YLL savings`/ econdata$irate
    econdata$`Absenteeism cost savings`<- econdata$`Absenteeism cost savings`/ econdata$irate
    econdata$`Presenteeism cost savings`<- econdata$`Presenteeism cost savings`/ econdata$irate
    econdata$irate <- NULL
    econdata$n <- NULL
    
    econdata<- econdata %>% summarise(across(everything(), list(sum)))
    econdata<- econdata/(10^240)
    econpie<-data.frame(c(econdata$`Screening cost`+econdata$`RFR cost`+econdata$`Diabetes treatment cost`+econdata$`Monitoring cost`,econdata$`YLL savings`+econdata$`Absenteeism cost savings`+econdata$`Presenteeism cost savings`))
    econpie$`Value`[1]<-"Costs"
    econpie$`Value`[2]<-"Benefits"
    colnames(econpie) <- c("Value","Variable")
    econdata$pop<-""
    econdata <- melt(econdata, id = "pop")
    econdata$pop<-NULL


  
    
    ###Create graphs  

    p3<-  ggplot() + geom_bar(data = piedata, aes(x = reorder(population, -value), y = value, fill = population), position = "dodge", stat = "identity") + theme(legend.position="none", axis.text.x = element_text(size = 10), axis.title.x = element_blank(),axis.title.y = element_blank()) +
      labs(title = "Population distribution at the end of the simulation (%)")
    
    p4<-  ggplot() + geom_bar(data = econdata, aes(x = variable, y = value, fill = variable), position = "dodge", stat = "identity") + theme(legend.position="none", axis.text.x = element_text(size = 8), axis.title.x = element_blank(),axis.title.y = element_blank()) + 
      labs(title = "NPV of selected economic outcomes (Billions)")

    
#    econpie <-econpie  %>% mutate(Variable = factor(Variable, levels = c("Costs","Benefits")), cumulative = cumsum(Value), midpoint = cumulative - Value / 2, label = paste0(Variable, " ", round(Value / sum(Value) * 100, 1), "%"))
 
    econpie %>%
      arrange(desc(Value)) %>%
      mutate(prop = percent(Value / sum(Value))) -> econpie 
    
    p5 <- ggplot(econpie, aes(x = "", y = Value, fill = fct_inorder(Variable))) + geom_bar(width = 1, stat = "identity") + coord_polar("y", start = 0) + geom_label_repel(aes(label = prop), size=5, show.legend = F, nudge_x = 2) + guides(fill = guide_legend(title = "Cost-benefit")) + labs(title = "Cost-benefit analysis (NPV)") + theme_void()
    
    
    
    #Arrangement of graphs
    grid.arrange(p3, p4, p5, layout_matrix = matrix(c(1, 2, 3, 3), nrow = 2))

  }) 
}

# Run the application
shinyApp(ui = ui, server = server)
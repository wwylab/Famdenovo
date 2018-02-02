require(LFSPRO)

Famdenovo <- function(family, cancer, person.id, mutation = NA, gene = "TP53") {
  
  if (gene != "TP53") {
    print("Famdenevoe can only analyze families with TP53 mutation.")
    return()
  }
  
  fam <- family[c("id", "fid", "mid", "gender", "age")]
  fam$fam.id <- 1
  can <- cancer[c("id", "cancer.type", "diag.age")]
  can$fam.id <- 1
  
  # family data: add indexes for each family and person
  fam1 <- fam
  count <- 0
  for (i in 1:length(unique(fam$fam.id))) {
    for (j in 1:sum(fam$fam.id == unique(fam$fam.id)[i])){
      count = count + 1
      fam1$index.fam[count]  <- i
      fam1$index.per[count]  <- j
    }
  }
  
  fam.tp <- merge(fam1, data.frame(fam.id=fam1$fam.id, fid=fam1$id, index.f=fam1$index.per), c("fam.id", "fid"), all.x=T)
  fam2 <- merge(fam.tp, data.frame(fam.id=fam1$fam.id, mid=fam1$id, index.m=fam1$index.per), c("fam.id", "mid"), all.x=T)
  fam2[is.na(fam2)] <- 0
  
  # cancer data: add indexes for family, person, and cancer type
  can1 <- merge(can, fam1[c("fam.id", "id", "index.fam", "index.per")], c("fam.id", "id"), all.x=T)
  LFSpro.cancer.type2 <- data.frame(cancer.type=names(LFSpro.cancer.type), index.can=LFSpro.cancer.type, stringsAsFactors=F)
  can2 <- merge(can1, LFSpro.cancer.type2, "cancer.type", all.x=T)
  
  # combine family and cancer data
  fam3 <- fam2[order(fam2$index.fam, fam2$index.per),]
  fam3 <- data.frame(fam.id=fam3$fam.id, id=fam3$index.per, fid=fam3$index.f, mid=fam3$index.m, gender=fam3$gender, age=fam3$age)
  can3 <- can2[order(can2$index.fam, can2$index.per),]
  can3 <- data.frame(fam.id=can3$fam.id, id=can3$index.per, cancer.type=can3$index.can, diag.age=can3$diag.age)
  fam.can <- CombineData(fam3, can3) 

  if (sum(is.na(mutation))) {
    print ("Warning: Famdenovo output is only applicatble to mutation carriers")
    counselees <- data.frame(fam.id=1, id=person.id, stringsAsFactors=F)
    counselees <- merge(counselees, fam1[c("fam.id", "id", "index.fam", "index.per")], c("fam.id", "id"), all.x=T)
  } else {
    mut <- mutation[c("id", "mut.state")]
    mut$fam.id <- 1
    
    # family and mutation data: add mutation state of father's and mother's
    fam.mut1 <- merge(fam2, data.frame(fam.id=mut$fam.id, id=mut$id, mut=mut$mut.state), c("fam.id", "id"), all.x=T)
    fam.mut.tp <- merge(fam.mut1, data.frame(fam.id=fam.mut1$fam.id, fid=fam.mut1$id, mut.f=fam.mut1$mut), c("fam.id", "fid"), all.x=T)
    fam.mut2 <- merge(fam.mut.tp, data.frame(fam.id=fam.mut1$fam.id, mid=fam.mut1$id, mut.m=fam.mut1$mut), c("fam.id", "mid"), all.x=T)
    
    # counselees: seperate by familial, de novo and unknown states and invalid person id
    fam.mut3 <- fam.mut2
    for (i in 1:nrow(fam.mut2)) {
      if(!is.na(fam.mut2$mut[i])) {
        if (fam.mut2$mut[i] %in% c("M")) {
          if (fam.mut2$mut.f[i] %in% c("M") || fam.mut2$mut.m[i] %in% c("M")) fam.mut3$state[i] <- "familial"
          else if(fam.mut2$mut.f[i] %in% c("W") && fam.mut2$mut.f[i] %in% c("W")) fam.mut3$state[i] <- "de novo"
          else fam.mut3$state[i] <- "unknown"
        }
        else
          fam.mut3[i, "state"] <- NA
      }
    }
    counselees <- fam.mut3[!is.na(fam.mut3$state), c("fam.id", "id", "index.fam", "index.per", "state")]
    person.id.invalid <- person.id[! person.id %in% counselees$id]
    if (length(person.id.invalid))
      print(paste("The following ids are not carriers:", paste(person.id.invalid, collapse = ", ")))
    counselees <- counselees[counselees$id %in% person.id,]
  }
  
  # calculate probability
  MAF <- 0.0003
  mRate <- 0.00006
  allef <- list(c(1 - MAF, MAF))
  
  pval <- numeric()
  for(i in 1:nrow(counselees)) {
    FaDFta <- fam.can[[counselees$index.fam[i]]]
    counselee.id <- counselees$index.per[i]
    ped <- data.frame(ID=FaDFta$id, Gender=FaDFta$gender, FatherID=FaDFta$fid, MotherID=FaDFta$mid)
    
    # modify likelihood
    lik <- calLK(FaDFta, LFSpenet.2010)
    lik_mod <- lik
    lik_mod[ped$ID == ped$FatherID[ped$ID == counselee.id], 2:3] <- 0
    lik_mod[ped$ID == ped$MotherID[ped$ID == counselee.id], 2:3] <- 0
    
    # P(Gc=1|Gm=0,Gf=0,D,P) = p_mod_counselee[2] + p_mod_counselee[3]
    p_mod_counselee <- peelingRC(allef, lik_mod, ped, counselee.id, 1, mRate)
    p_mod_Gc1 <- p_mod_counselee[2] + p_mod_counselee[3]
    
    # P(Gc=1|D,P)
    p_counselee <- peelingRC(allef, lik, ped, counselee.id, 1, mRate)
    p_Gc1 <- p_counselee[2] + p_counselee[3]
    
    # P(Gf=0|D,P)
    tp_id <- ped$FatherID[ped$ID == counselee.id]
    p_father <- peelingRC(allef, lik, ped, tp_id, 1, mRate)
    p_Gf0 <- p_father[1]
    
    # P(Gm=0|D,P)
    tp_id <- ped$MotherID[ped$ID == counselee.id]
    p_mother <- peelingRC(allef, lik, ped, tp_id, 1, mRate)
    p_Gm0 <- p_mother[1]
    
    # P(Gc is denovo|Gc is a germline, D, P) = P(Gm=0|D,P) * P(Gf=0|D,P) * P(Gc=1|Gm=0,Gf=0,D,P) / P(Gc=1|D,P)
    pval[i] <-  p_Gm0 * p_Gf0 * p_mod_Gc1/ p_Gc1
  }
  
  pred_prob <- data.frame(id=counselees$id, prob.denovo=pval)
  return(pred_prob)
}


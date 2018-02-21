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
  count <- 1
  fams <- unique(fam$fam.id)
  for (i in 1:length(fams)) {
    for (j in 1:sum(fam$fam.id == fams[i])){
      fam1$index.fam[count]  <- i
      fam1$index.per[count]  <- j
      count = count + 1
    }
  }
  
  # family data: add indexes for father and mother
  fam2 <- merge(fam1, data.frame(fam.id=fam1$fam.id, fid=fam1$id, index.f=fam1$index.per), c("fam.id", "fid"), all.x=T)
  fam2 <- merge(fam2, data.frame(fam.id=fam1$fam.id, mid=fam1$id, index.m=fam1$index.per), c("fam.id", "mid"), all.x=T)
  fam2[is.na(fam2)] <- 0
  
  # cancer data: add indexes for family, person, and cancer type
  can1 <- merge(can, fam1[c("fam.id", "id", "index.fam", "index.per")], c("fam.id", "id"), all.x=T)
  can2 <- merge(can1, LFSpro.cancer.type, "cancer.type", all.x=T)
  
  # combine family and cancer data
  fam3 <- data.frame(fam.id=fam2$fam.id, id=fam2$index.per, fid=fam2$index.f, mid=fam2$index.m, gender=fam2$gender, age=fam2$age, stringsAsFactors=F)
  can3 <- data.frame(fam.id=can2$fam.id, id=can2$index.per, cancer.type=can2$index.can, diag.age=can2$diag.age)
  fam.can <- CombineData(fam3, can3)

  # mutation data
  if (sum(is.na(mutation))) {
    print ("Warning: Famdenovo output is only applicatble to mutation carriers")
    counselees <- data.frame(fam.id=1, id=person.id, stringsAsFactors=F)
    counselees <- merge(counselees, fam1[c("fam.id", "id", "index.fam", "index.per")], c("fam.id", "id"), all.x=T)
  } else {
    mut <- mutation[c("id", "mut.state")]
    mut$fam.id <- 1
    
    # family and mutation data: add mutation states of father and mother
    fam.mut1 <- merge(fam2, data.frame(fam.id=mut$fam.id, id=mut$id, mut=mut$mut.state), c("fam.id", "id"), all.x=T)
    fam.mut2 <- merge(fam.mut1, data.frame(fam.id=fam.mut1$fam.id, fid=fam.mut1$id, mut.f=fam.mut1$mut), c("fam.id", "fid"), all.x=T)
    fam.mut2 <- merge(fam.mut2, data.frame(fam.id=fam.mut1$fam.id, mid=fam.mut1$id, mut.m=fam.mut1$mut), c("fam.id", "mid"), all.x=T)
    
    # counselees: denovo status
    fam.mut3 <- fam.mut2
    for (i in 1:nrow(fam.mut2)) {
      if (fam.mut2$mut[i] %in% "M") {
        if (fam.mut2$mut.f[i] %in% "M" || fam.mut2$mut.m[i] %in% "M") {
          fam.mut3$state[i] <- "familial"
        } else if (fam.mut2$mut.f[i] %in% "W" && fam.mut2$mut.m[i] %in% "W") {
          fam.mut3$state[i] <- "denovo"
        } else {
          fam.mut3$state[i] <- "unknown"
        }
      } else {
        fam.mut3$state[i] <- NA
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
    counselee.id <- counselees$index.per[i]
    for(j in 1:length(fam.can)) {
      if (unique(fam.can[[j]]$fam.id) == counselees$fam.id[i]) break
    }
    FaDFta <- fam.can[[j]]
    ped <- data.frame(ID=FaDFta$id, Gender=FaDFta$gender, FatherID=FaDFta$fid, MotherID=FaDFta$mid)
    
    # modify likelihood
    lik <- calLK(FaDFta, LFSpenet.2010)
    lik.mod <- lik
    lik.mod[ped$ID == ped$FatherID[ped$ID == counselee.id], 2:3] <- 0
    lik.mod[ped$ID == ped$MotherID[ped$ID == counselee.id], 2:3] <- 0
    
    # P(Gc=1|D,P)
    p.counselee <- peelingRC(allef, lik, ped, counselee.id, 1, mRate)
    p.Gc1 <- p.counselee[2] + p.counselee[3]
    
    # P(Gc=1|Gm=0,Gf=0,D,P)
    p.counselee.mod <- peelingRC(allef, lik.mod, ped, counselee.id, 1, mRate)
    p.Gc1.mod <- p.counselee.mod[2] + p.counselee.mod[3]
    
    # P(Gf=0|D,P)
    fid <- ped$FatherID[ped$ID == counselee.id]
    if (fid == 0) {
      p.Gf0 <- 1- (1 - MAF)^2
    } else {
      p.father <- peelingRC(allef, lik, ped, fid, 1, mRate)
      p.Gf0 <- p.father[1]
    }
    
    # P(Gm=0|D,P)
    mid <- ped$MotherID[ped$ID == counselee.id]
    if (mid == 0) {
      p.Gm0 <- 1- (1 - MAF)^2
    } else {
      p.mother <- peelingRC(allef, lik, ped, mid, 1, mRate)
      p.Gm0 <- p.mother[1]
    }
    
    # P(Gc is denovo|Gc is a germline, D, P) = P(Gm=0|D,P) * P(Gf=0|D,P) * P(Gc=1|Gm=0,Gf=0,D,P) / P(Gc=1|D,P)
    pval[i] <-  p.Gm0 * p.Gf0 * p.Gc1.mod/ p.Gc1
  }
  
  pred.prob <- data.frame(id=counselees$id, prob.denovo=pval); pred.prob
  return(pred.prob)
}


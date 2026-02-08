############################################################
##  MICOM-style co-operative trade-off for BacArena/Sybil ##
############################################################
#
#  build_comm_model()     – combine all unique species models
#  ctFBA()                – community FBA
#  apply_ctFBA_to_arena() – writes mu_i,target back into the arena
#                           by capping each species’ biomass flux
# ---------------------------------------------------------

library(sybil)
library(Matrix)
library(glpkAPI)
library(sybilSBML) 
library(BacArena)

# ------------------------------------------------------------------------------------------
# sync_metabolite_slots() - helper to combine multiple models into single 'community' model
# ------------------------------------------------------------------------------------------
sync_metabolite_slots <- function(mod) {
  
  S_matrix <- S(mod)
  n_mets <- nrow(S_matrix)
  current_rownames <- rownames(S_matrix)
  
  cat("Syncing metabolite slots: S matrix has", n_mets, "rows\n")
  
  # Handle met_id slot - make it exactly match the S matrix rownames
  if (!is.null(current_rownames) && length(current_rownames) == n_mets) {
    # Use the rownames directly
    met_id(mod) <- current_rownames
  } else {
    # Generate generic names if rownames are missing/wrong
    generic_names <- paste0("met_", seq_len(n_mets))
    met_id(mod) <- generic_names
    if (!is.null(current_rownames)) {
      rownames(S(mod)) <- generic_names
    }
  }
  
  # Handle met_name slot - make it same length as met_id
  current_met_names <- met_name(mod)
  target_length <- length(met_id(mod))
  
  if (length(current_met_names) != target_length) {
    # Extend or truncate to match
    if (length(current_met_names) < target_length) {
      # Extend with empty strings or copy from met_id
      extension <- met_id(mod)[(length(current_met_names) + 1):target_length]
      met_name(mod) <- c(current_met_names, extension)
    } else {
      # Truncate
      met_name(mod) <- current_met_names[1:target_length]
    }
  }
  
  # Handle met_comp slot if it exists
  if ("met_comp" %in% slotNames(mod)) {
    current_met_comp <- met_comp(mod)
    if (length(current_met_comp) != target_length) {
      if (length(current_met_comp) < target_length) {
        # Extend with compartment 1
        extension <- rep(1L, target_length - length(current_met_comp))
        met_comp(mod) <- c(current_met_comp, extension)
      } else {
        # Truncate
        met_comp(mod) <- current_met_comp[1:target_length]
      }
    }
  }
  
  # Verify the final state
  final_met_id_length <- length(met_id(mod))
  final_met_name_length <- length(met_name(mod))
  
  cat("Final lengths - S matrix:", n_mets, 
      "met_id:", final_met_id_length, 
      "met_name:", final_met_name_length, "\n")
  
  if (final_met_id_length != n_mets || final_met_name_length != n_mets) {
    stop("Failed to sync metabolite slots properly!")
  }
  
  return(mod)
}


# -----------------------------------------------------------------
# Helper function to ensure combined S matrix has proper rownames
# -----------------------------------------------------------------
fix_matrix_rownames <- function(mod) {
  S_matrix <- S(mod)
  
  # Check if rownames are missing
  if (is.null(rownames(S_matrix))) {
    # Get metabolite IDs from the model
    met_ids <- met_id(mod)
    
    # Ensure we have the right number of metabolite IDs
    if (length(met_ids) == nrow(S_matrix)) {
      rownames(S_matrix) <- met_ids
      S(mod) <- S_matrix
      # cat("Fixed rownames for model with", nrow(S_matrix), "metabolites\n")
    } else {
      # Generate generic metabolite names if met_id doesn't match
      generic_names <- paste0("met_", 1:nrow(S_matrix))
      rownames(S_matrix) <- generic_names
      S(mod) <- S_matrix
      
      # Also update met_id slot to match
      met_id(mod) <- generic_names
      cat("Generated", length(generic_names), "generic metabolite names\n")
    }
  }
  
  # Also ensure colnames exist for reactions
  if (is.null(colnames(S_matrix))) {
    react_ids <- react_id(mod)
    if (length(react_ids) == ncol(S_matrix)) {
      colnames(S(mod)) <- react_ids
    } else {
      colnames(S(mod)) <- paste0("rxn_", 1:ncol(S_matrix))
    }
  }
  
  return(mod)
}




# -----------------------------------------------------------------
# build_comm_model() - handle Matrix types properly
# -----------------------------------------------------------------
build_comm_model <- function(arena, abund_weight = TRUE) {
  
  # Collect unique model objects
  specs <- arena@specs
  uniq  <- !duplicated(vapply(specs, function(x) x@model@mod_desc, ""))
  mods  <- lapply(specs[uniq], slot, "model")
  
  cat("Building community model from", length(mods), "unique species\n")
  
  # Fix rownames for all models
  cat("Fixing matrix rownames...\n")
  mods <- lapply(mods, fix_matrix_rownames)
  
  # make reaction IDs unique before merging models
  mods <- lapply(seq_along(mods), function(i) {
    tag <- paste0("sp", i)                    # or use mod@mod_desc
    unique_react_ids(mods[[i]], tag)
  })
  
  # Start with first model as base
  comm <- mods[[1]]
  
  # Ensure the base model has a numeric sparse matrix
  S_base <- S(comm)
  if (class(S_base)[1] == "ngCMatrix") {
    S(comm) <- as(S_base, "dgCMatrix")
    cat("Converted base model to numeric sparse matrix\n")
  }
  
  if (length(mods) > 1) {
    for (k in 2:length(mods)) {
      m <- mods[[k]]
      
      # Get metabolite names
      comm_mets <- rownames(S(comm))
      m_mets <- rownames(S(m))
      all_mets <- union(comm_mets, m_mets)
      
      cat("Combining model", k, "- total metabolites:", length(all_mets), "\n")
      
      # Convert to regular matrices to avoid Matrix type issues
      S_comm <- as.matrix(S(comm))
      S_m <- as.matrix(S(m))
      
      # Pad comm matrix
      missing_in_comm <- setdiff(all_mets, comm_mets)
      if (length(missing_in_comm) > 0) {
        zero_rows_comm <- matrix(0.0, nrow = length(missing_in_comm), ncol = ncol(S_comm))
        rownames(zero_rows_comm) <- missing_in_comm
        # Ensure column names match
        if (!is.null(colnames(S_comm))) {
          colnames(zero_rows_comm) <- colnames(S_comm)
        }
        S_comm <- rbind(S_comm, zero_rows_comm)
      }
      
      # Pad m matrix  
      missing_in_m <- setdiff(all_mets, m_mets)
      if (length(missing_in_m) > 0) {
        zero_rows_m <- matrix(0.0, nrow = length(missing_in_m), ncol = ncol(S_m))
        rownames(zero_rows_m) <- missing_in_m
        # Ensure column names match
        if (!is.null(colnames(S_m))) {
          colnames(zero_rows_m) <- colnames(S_m)
        }
        S_m <- rbind(S_m, zero_rows_m)
      }
      
      # Reorder both to same metabolite order
      S_comm <- S_comm[all_mets, , drop = FALSE]
      S_m <- S_m[all_mets, , drop = FALSE]
      
      # Combine horizontally - ensure no column name conflicts
      combined_colnames <- c(colnames(S_comm), colnames(S_m))
      S_combined <- cbind(S_comm, S_m)
      colnames(S_combined) <- combined_colnames
      
      # Convert back to numeric sparse matrix explicitly
      S(comm) <- Matrix::Matrix(S_combined, sparse = TRUE)
      S(comm) <- as(S(comm), "dgCMatrix")  # Ensure numeric sparse matrix
      
      # Update reaction slots
      react_id(comm)  <- c(react_id(comm),  react_id(m))
      lowbnd(comm)    <- c(lowbnd(comm),    lowbnd(m))
      uppbnd(comm)    <- c(uppbnd(comm),    uppbnd(m))
      obj_coef(comm)  <- c(obj_coef(comm),  obj_coef(m))
      react_rev(comm) <- c(react_rev(comm), react_rev(m))
    }
  }
  
  slot(comm, "react_num") <- as.integer(ncol(S(comm)))
  slot(comm, "met_num")   <- as.integer(nrow(S(comm))) 
  
  nRxn <- ncol(S(comm))
  nMet <- nrow(S(comm))
  
  
  # COMPREHENSIVE SLOT SYNCHRONIZATION:
  align_cols <- function(A, B, fill = NA) {
    if (is.null(colnames(A))) colnames(A) <- seq_len(ncol(A))
    if (is.null(colnames(B))) colnames(B) <- seq_len(ncol(B))
    
    missA <- setdiff(colnames(B), colnames(A))
    missB <- setdiff(colnames(A), colnames(B))
    
    if (length(missA))
      A <- cbind(A,
                 matrix(fill, nrow=nrow(A), ncol=length(missA),
                        dimnames = list(NULL, missA)))
    if (length(missB))
      B <- cbind(B,
                 matrix(fill, nrow=nrow(B), ncol=length(missB),
                        dimnames = list(NULL, missB)))
    
    A <- A[, sort(colnames(A)), drop = FALSE]
    B <- B[, sort(colnames(B)), drop = FALSE]
    list(A=A, B=B)
  }
  
  # ===== VECTOR SLOTS =====
  all_reaction_slots <- c("react_id","react_name","react_rev","lowbnd","uppbnd","obj_coef",
                          "react_single","react_de","gprRules","gpr")
  all_metabolite_slots <- c("met_id","met_name","met_comp","met_single","met_de")
  
  for (sn in all_reaction_slots) {
    if (sn %in% slotNames(comm)) {
      current_length <- length(slot(comm, sn))
      if (current_length != nRxn) {
        length(slot(comm, sn)) <- nRxn
      }
    }
  }
  
  for (sn in all_metabolite_slots) {
    if (sn %in% slotNames(comm)) {
      current_length <- length(slot(comm, sn))
      if (current_length != nMet) {
        length(slot(comm, sn)) <- nMet
      }
    }
  }
  
  # ===== MATRIX SLOTS =====
  # Handle react_attr - ensure it has nRxn rows
  if ("react_attr" %in% slotNames(comm)) {
    current_attr <- react_attr(comm)
    if (nrow(current_attr) != nRxn) {
      if (nrow(current_attr) < nRxn) {
        # Add missing rows with NA values, ensuring column names match
        missing_rows <- nRxn - nrow(current_attr)
        if (ncol(current_attr) > 0) {
          new_rows <- matrix(NA, nrow = missing_rows, ncol = ncol(current_attr))
          colnames(new_rows) <- colnames(current_attr)
          react_attr(comm) <- rbind(current_attr, new_rows)
        }
      } else {
        # Truncate to correct size
        react_attr(comm) <- current_attr[1:nRxn, , drop = FALSE]
      }
    }
  }
  
  # Handle met_attr - ensure it has nMet rows  
  if ("met_attr" %in% slotNames(comm)) {
    current_attr <- met_attr(comm)
    if (nrow(current_attr) != nMet) {
      if (nrow(current_attr) < nMet) {
        # Add missing rows with NA values, ensuring column names match
        missing_rows <- nMet - nrow(current_attr)
        if (ncol(current_attr) > 0) {
          new_rows <- matrix(NA, nrow = missing_rows, ncol = ncol(current_attr))
          colnames(new_rows) <- colnames(current_attr)
          met_attr(comm) <- rbind(current_attr, new_rows)
        }
      } else {
        # Truncate to correct size
        met_attr(comm) <- current_attr[1:nMet, , drop = FALSE]
      }
    }
  }
  
  
  # Handle subSys - ensure proper dimensions
  if ("subSys" %in% slotNames(comm)) {
    current_subsys <- subSys(comm)
    if (is(current_subsys, "Matrix")) {
      # Check if reactions are rows or columns and adjust accordingly
      if (nrow(current_subsys) == nRxn || abs(nrow(current_subsys) - nRxn) <= 1) {
        # Reactions are rows - adjust rows to match nRxn
        if (nrow(current_subsys) < nRxn) {
          new_rows <- Matrix::sparseMatrix(i = integer(0), j = integer(0),
                                           dims = c(nRxn - nrow(current_subsys), ncol(current_subsys)))
          subSys(comm) <- rbind(current_subsys, new_rows)
        } else if (nrow(current_subsys) > nRxn) {
          subSys(comm) <- current_subsys[1:nRxn, , drop = FALSE]
        }
      } else if (ncol(current_subsys) == nRxn || abs(ncol(current_subsys) - nRxn) <= 1) {
        # Reactions are columns - adjust columns to match nRxn
        if (ncol(current_subsys) < nRxn) {
          new_cols <- Matrix::sparseMatrix(i = integer(0), j = integer(0),
                                           dims = c(nrow(current_subsys), nRxn - ncol(current_subsys)))
          subSys(comm) <- cbind(current_subsys, new_cols)
        } else if (ncol(current_subsys) > nRxn) {
          subSys(comm) <- current_subsys[, 1:nRxn, drop = FALSE]
        }
      } else {
        # Create new empty subSys matrix
        cat("subSys dimensions don't match - creating new matrix\n")
        subSys(comm) <- Matrix::sparseMatrix(i = integer(0), j = integer(0),
                                             dims = c(nRxn, 0))
      }
    } else {
      # Create new subSys matrix if it's not a Matrix object
      subSys(comm) <- Matrix::sparseMatrix(i = integer(0), j = integer(0),
                                           dims = c(nRxn, 0))
    }
  }
  
  # Update model description
  mod_desc(comm) <- paste("COMM", length(mods), "species")
  
  # Sync metabolite slots
  comm <- sync_metabolite_slots(comm)
  
  # Find biomass reactions
  bio_pos <- c()
  bio_ids <- c()
  current_offset <- 0
  
  for (i in 1:length(mods)) {
    m <- mods[[i]]
    bio_candidates <- grep("bio|biomass", react_id(m), ignore.case = TRUE)
    
    if (length(bio_candidates) > 0) {
      bio_local <- bio_candidates[1]
      bio_global <- current_offset + bio_local
      bio_pos <- c(bio_pos, bio_global)
      bio_ids <- c(bio_ids, react_id(m)[bio_local])
    } else {
      warning("No biomass reaction found for model ", i)
    }
    
    current_offset <- current_offset + react_num(m)
  }
  
  cat("Found", length(bio_pos), "biomass reactions\n")
  
  # ===== CRITICAL FIX: Handle gene-related slots properly =====
  
  # Collect all genes from original models to maintain proper structure
  all_genes <- character(0)
  all_gpr_rules <- character(0)
  all_rxn_gene_mats <- list()
  
  # Collect gene information from all models
  for (i in 1:length(mods)) {
    m <- mods[[i]]
    
    if ("genes" %in% slotNames(m)) {
      model_genes <- genes(m)
      if (is.list(model_genes)) {
        # Extract unique genes from all reactions in this model
        model_gene_chars <- unique(unlist(model_genes))
        all_genes <- unique(c(all_genes, model_gene_chars))
      }
    }
    
    if ("gprRules" %in% slotNames(m)) {
      all_gpr_rules <- c(all_gpr_rules, gprRules(m))
    }
    
    if ("rxnGeneMat" %in% slotNames(m)) {
      all_rxn_gene_mats[[i]] <- rxnGeneMat(m)
    }
  }
  
  # Set gene_num based on collected genes
  total_genes <- length(all_genes)
  if ("gene_num" %in% slotNames(comm)) {
    slot(comm, "gene_num") <- as.integer(total_genes)
  }
  
  # Create proper genes slot structure
  if ("genes" %in% slotNames(comm)) {
    # Create a list of length nRxn, each with empty character vector
    genes_list <- vector("list", nRxn)
    for (i in 1:nRxn) {
      genes_list[[i]] <- character(0)
    }
    genes(comm) <- genes_list
  }
  
  # Set gprRules 
  if ("gprRules" %in% slotNames(comm)) {
    if (length(all_gpr_rules) == nRxn) {
      gprRules(comm) <- all_gpr_rules
    } else {
      gprRules(comm) <- rep("", nRxn)
    }
  }
  
  # Set gpr slot
  if ("gpr" %in% slotNames(comm)) {
    gpr(comm) <- rep("", nRxn)
  }
  
  # Create rxnGeneMat with proper dimensions: nRxn x total_genes
  if ("rxnGeneMat" %in% slotNames(comm)) {
    rxnGeneMat(comm) <- Matrix::sparseMatrix(i = integer(0), j = integer(0),
                                             dims = c(nRxn, total_genes))
  }
  
  # Set allGenes slot if it exists
  if ("allGenes" %in% slotNames(comm)) {
    allGenes(comm) <- all_genes
  }
  
  # Handle subSys if not already handled above
  if ("subSys" %in% slotNames(comm)) {
    if (!is(subSys(comm), "Matrix") || 
        (nrow(subSys(comm)) != nRxn && ncol(subSys(comm)) != nRxn)) {
      subSys(comm) <- Matrix::sparseMatrix(i = integer(0), j = integer(0),
                                           dims = c(nRxn, 0))
    }
  }
  
  # Final validation
  tryCatch({
    validObject(comm)
    cat("Community model validation successful\n")
  }, error = function(e) {
    stop("Community model validation failed: ", e$message)
  })
  
  return(list(model = comm, bio_pos = bio_pos, bio_ids = bio_ids))
}




### =============================================================
# Main ctFBA function that works with community models
### =============================================================
##  comm_result         : modelorg object that already contains
##                 – individual biomass reactions
##
##  abund_vec     : numeric vector of *initial abundances* (same order as bio_pos)
##  tau           : trade-off (0–1) — keep ≥ tau*mu_Cmax
##  min_growth    : minimal mu_i cap if analytic formula gives < min_growth
##  solver        : "glpkAPI" (default) or "cplexAPI"
##  verbose       : logical
##

ctFBA <- function(comm_result, 
                  tau = 0.5, 
                  abund_vec = NULL,
                  min_growth = 1e-6,
                  solver = "glpkAPI",
                  verbose = TRUE) {
  
  # Extract components from community model result
  if (is.list(comm_result) && "model" %in% names(comm_result)) {
    comm <- comm_result$model
    bio_pos <- comm_result$bio_pos
    
    # If no abundance vector provided, use equal weights
    if (is.null(abund_vec)) {
      abund_vec <- rep(1, length(bio_pos))
    }
  } else {
    # Assume comm_result is the model directly
    comm <- comm_result
    bio_pos <- grep("bio|biomass", react_id(comm), ignore.case = TRUE)
    # bio_pos <- bio_pos[!grepl("COMM_BIOMASS", react_id(comm)[bio_pos])]  # exclude community biomass
    
    if (is.null(abund_vec)) {
      abund_vec <- rep(1, length(bio_pos))
    }
  }
  
  if (length(bio_pos) != length(abund_vec)) {
    stop("bio_pos and abund_vec must have the same length")
  }
  
  if (verbose) {
    cat("ctFBA: tau =", tau, ", biomass reactions =", length(bio_pos), "\n")
  }
  
  ## helper function just in case: sanitize model for LP solver ----------------------
  sanitize_for_lp <- function(mod, lb_def = 0, ub_def = 1000) {
    lb <- lowbnd(mod); lb[is.na(lb)] <- lb_def; lowbnd(mod) <- lb
    ub <- uppbnd(mod); ub[is.na(ub)] <- ub_def; uppbnd(mod) <- ub
    rv <- react_rev(mod); rv[is.na(rv)] <- FALSE; react_rev(mod) <- rv
    oc <- obj_coef(mod); oc[is.na(oc)] <- 0; obj_coef(mod) <- oc
    
    # Fix inconsistent bounds
    inconsistent <- lb > ub
    if (any(inconsistent)) {
      ub[inconsistent] <- lb[inconsistent]
      uppbnd(mod) <- ub
    }
    
    mod
  }
  
  comm <- sanitize_for_lp(comm)
  
  ## ---------- STEP 0 : set weighted objective -----------------
  rel_abund <- abund_vec / sum(abund_vec)          # a[i] 
  
  obj_coef(comm) <- rep(0, react_num(comm))        # reset everything
  obj_coef(comm)[bio_pos] <- rel_abund             # weight each biomass rxn
  
  ## ---------- STEP 1: Maximize community biomass ----------------------
  # Reset objectives and set COMM_BIOMASS = 1
  # obj_coef(comm) <- rep(0, length(react_id(comm)))
  # comm <- changeObjFunc(comm, react = "COMM_BIOMASS", obj_coef = 1)
  
  # Solve
  tryCatch({
    sol1 <- sybil::optimizeProb(comm, solver = solver)
    
    if (lp_ok(sol1) != 0 || lp_stat(sol1) != 5) {
      if (verbose) cat("Warning: LP failed when maximizing community biomass\n")
      return(rep(min_growth, length(bio_pos)))  # fallback
    }
    
    muCmax <- lp_obj(sol1) # now a weighted mean
    if (verbose) cat("µCmax =", muCmax, "\n")
    
  }, error = function(e) {
    if (verbose) cat("Error in community biomass optimization:", e$message, "\n")
    return(rep(min_growth, length(bio_pos)))  # fallback
  })
  
  ## ---------- STEP 2: Analytic L-2 balanced growth rates --------------
  rel_abund <- abund_vec / sum(abund_vec)             # = aᵢ (relative abundances) normalized weights
  scale <- (tau * muCmax) / sum(rel_abund^2)          # = (τ * μc_max) / Σᵢ aᵢ²  - scaling factor
  mu_star <- pmax(rel_abund * scale, min_growth)      # = aᵢ * (τ * μc_max) / Σᵢ aᵢ² - balanced growth rates
  
  if (verbose) {
    cat("Balanced growth rates: min =", min(mu_star), ", max =", max(mu_star), "\n")
  }
  
  # Return just the growth rates (for compatibility with choose_tau)
  return(mu_star)
}


apply_ctFBA_to_arena <- function(arena,
                                 tau      = 0.5,
                                 solver   = "glpkAPI",
                                 verbose  = TRUE,
                                 type_specie_table,
                                 min_growth = 1e-6) {
  
  if (verbose) cat("Applying ctFBA to arena with tau =", tau, "\n")
  
  ## -------- 1. build community & compute caps -----------------
  comm_result <- build_comm_model(arena)
  
  ## current absolute counts per *Type*
  type_counts <- table(arena@orgdat$type)         # names = "1","2",…
  spec_names  <- names(arena@specs)               # model descriptors
  type_ids    <- type_specie_table$Type[
    match(spec_names, type_specie_table$Name)]
  
  abund_vec <- as.numeric(type_counts[as.character(type_ids)])
  abund_vec[is.na(abund_vec)] <- 0                # species absent this hour
  abund_vec <- abund_vec / sum(abund_vec)         # relative abundance
  
  mu_star <- ctFBA(comm_result,
                   tau        = tau,
                   abund_vec  = abund_vec,
                   solver     = solver,
                   verbose    = verbose)
  
  ## ---- name the caps by Type so we can look them up ----------
  names(mu_star) <- as.character(type_ids)        # e.g. "15"  "8" …
  
  ## -------- 2. apply caps to every species now present --------
  for (sp in seq_along(arena@specs)) {
    
    sp_name <- names(arena@specs)[sp]           # descriptor
    type_id <- type_specie_table$Type[
      match(sp_name, type_specie_table$Name)]
    
    ## cap for this species (or fallback)
    mu_cap  <- mu_star[as.character(type_id)]
    if (is.na(mu_cap) || mu_cap <= 0)
      mu_cap <- min_growth                    # newcomers / extinct
    
    ## set UB of *all* biomass reactions in its GEM
    mod_k   <- arena@specs[[sp]]@model
    bio_idx <- grep("bio|biomass", react_id(mod_k), ignore.case = TRUE)
    
    if (length(bio_idx)) {
      upp            <- uppbnd(mod_k)
      upp[bio_idx]   <- mu_cap + 1e-8         # small slack
      uppbnd(mod_k)  <- upp
      arena@specs[[sp]]@model <- mod_k
    }
  }
  
  if (verbose)
    cat("Applied growth caps:\n",
        paste(sprintf("%s:%.3g", names(mu_star), mu_star), collapse = "  "),
        "\n")
  
  return(arena)
}


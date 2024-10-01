% This code was created by Loïc Marrec 

function pfix = Moran_dioecy(Nrep, delta, hF, hM, NAAF0, NABF0, NBBF0, NAAM0, NABM0, NBBM0)
                                                   
        pfix = 0;
        
        % Loop over the number of replicates
        
        for i = 1 : Nrep

            % Initialization

            NAAF = NAAF0;   % Number of AA females          
            NABF = NABF0;   % Number of AB females
            NBBF = NBBF0;   % Number of BB females
            NAAM = NAAM0;   % Number of AA males          
            NABM = NABM0;   % Number of AB males
            NBBM = NBBM0;   % Number of BB males
            
            % The simulation continues until only one genotype remains in the population

            while not((NAAF+NAAM ~= 0 && NABF+NABM == 0 && NBBF+NBBM == 0) || (NAAF+NAAM == 0 && NABF+NABM == 0 && NBBF+NBBM ~= 0))
                
                % Compute the probabilities of picking one AA, AB, or BB
                % individidual
                
                bAAF = (1+delta)*NAAF/((1+delta)*NAAF+(1+hF*delta)*NABF+NBBF);
                bABF = (1+hF*delta)*NABF/((1+delta)*NAAF+(1+hF*delta)*NABF+NBBF);
                bBBF = NBBF/((1+delta)*NAAF+(1+hF*delta)*NABF+NBBF);
                bAAM = (1+delta)*NAAM/((1+delta)*NAAM+(1+hM*delta)*NABM+NBBM);
                bABM = (1+hM*delta)*NABM/((1+delta)*NAAM+(1+hM*delta)*NABM+NBBM);
                bBBM = NBBM/((1+delta)*NAAM+(1+hM*delta)*NABM+NBBM);
                
                % Compute the probabilities of picking one A or B gamete
                
                bAF = bAAF+bABF/2;
                bBF = bBBF+bABF/2;
                bAM = bAAM+bABM/2;
                bBM = bBBM+bABM/2;
                
                % Compute the probabilities that the offspring is AA, AB,
                % or AB
                
                bAA_off = bAF*bAM;
                bAB_off = bAF*bBM+bAM*bBF;
                bBB_off = bBF*bBM;

                % Generate a new offspring at random

                r1 = rand;
                r2 = rand;
                
                if r1 <= bAA_off/(bAA_off+bAB_off+bBB_off)
                    
                    if r2 <= .5
                    
                        % Pick randomly one individual for death
                
                        r3 = rand;
                        
                        if r3 <= NAAF/(NAAF+NABF+NBBF)

                            NAAF = NAAF-1;

                        elseif r3 > NAAF/(NAAF+NABF+NBBF) && r3 <= (NAAF+NABF)/(NAAF+NABF+NBBF)

                            NABF = NABF-1;

                        elseif r3 > (NAAF+NABF)/(NAAF+NABF+NBBF) && r3 <= 1

                            NBBF = NBBF-1;

                        end
                        
                        % A new offspring is born
                        
                        NAAF = NAAF+1;
                        
                    else
                        
                        % Pick randomly one individual for death
                
                        r3 = rand;
                        
                        if r3 <= NAAM/(NAAM+NABM+NBBM)

                            NAAM = NAAM-1;

                        elseif r3 > NAAM/(NAAM+NABM+NBBM) && r3 <= (NAAM+NABM)/(NAAM+NABM+NBBM)

                            NABM = NABM-1;

                        elseif r3 > (NAAM+NABM)/(NAAM+NABM+NBBM) && r3 <= 1

                            NBBM = NBBM-1;

                        end
                        
                        % A new offspring is born
                        
                        NAAM = NAAM+1;
                        
                    end
                    
                elseif r1 > bAA_off/(bAA_off+bAB_off+bBB_off) && r1 <= (bAA_off+bAB_off)/(bAA_off+bAB_off+bBB_off)
                    
                    if r2 <= .5
                    
                        % Pick randomly one individual for death
                
                        r3 = rand;
                        
                        if r3 <= NAAF/(NAAF+NABF+NBBF)

                            NAAF = NAAF-1;

                        elseif r3 > NAAF/(NAAF+NABF+NBBF) && r3 <= (NAAF+NABF)/(NAAF+NABF+NBBF)

                            NABF = NABF-1;

                        elseif r3 > (NAAF+NABF)/(NAAF+NABF+NBBF) && r3 <= 1

                            NBBF = NBBF-1;

                        end
                        
                        % A new offspring is born
                        
                        NABF = NABF+1;
                        
                    else
                        
                        % Pick randomly one individual for death
                
                        r3 = rand;
                        
                        if r3 <= NAAM/(NAAM+NABM+NBBM)

                            NAAM = NAAM-1;

                        elseif r3 > NAAM/(NAAM+NABM+NBBM) && r3 <= (NAAM+NABM)/(NAAM+NABM+NBBM)

                            NABM = NABM-1;

                        elseif r3 > (NAAM+NABM)/(NAAM+NABM+NBBM) && r3 <= 1

                            NBBM = NBBM-1;

                        end
                        
                        % A new offspring is born
                        
                        NABM = NABM+1;
                        
                    end
                    
                elseif r1 > (bAA_off+bAB_off)/(bAA_off+bAB_off+bBB_off) && r1 <= 1
                    
                    if r2 <= .5
                    
                        % Pick randomly one individual for death
                
                        r3 = rand;
                        
                        if r3 <= NAAF/(NAAF+NABF+NBBF)

                            NAAF = NAAF-1;

                        elseif r3 > NAAF/(NAAF+NABF+NBBF) && r3 <= (NAAF+NABF)/(NAAF+NABF+NBBF)

                            NABF = NABF-1;

                        elseif r3 > (NAAF+NABF)/(NAAF+NABF+NBBF) && r3 <= 1

                            NBBF = NBBF-1;

                        end
                        
                        % A new offspring is born
                        
                        NBBF = NBBF+1;
                        
                    else
                        
                        % Pick randomly one individual for death
                
                        r3 = rand;
                        
                        if r3 <= NAAM/(NAAM+NABM+NBBM)

                            NAAM = NAAM-1;

                        elseif r3 > NAAM/(NAAM+NABM+NBBM) && r3 <= (NAAM+NABM)/(NAAM+NABM+NBBM)

                            NABM = NABM-1;

                        elseif r3 > (NAAM+NABM)/(NAAM+NABM+NBBM) && r3 <= 1

                            NBBM = NBBM-1;

                        end
                        
                        % A new offspring is born
                        
                        NBBM = NBBM+1;
                        
                    end
                    
                end

            end
            
            % Check whether the mutation took over the population
            
            if NAAF+NAAM ~= 0
                
                pfix = pfix+1;
                
            end
        
        end
        
        % Calculate the fixation probability
        
        pfix = pfix/Nrep;
   
end
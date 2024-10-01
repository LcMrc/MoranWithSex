% This code was created by Loïc Marrec

function pfix = Moran_monoecy(Nrep, delta, h, NAA0, NAB0, NBB0)
                                                   
        pfix = 0;
        
        % Loop over the number of replicates
        
        for i = 1 : Nrep

            % Initialization

            NAA = NAA0;     % Number of AA individuals          
            NAB = NAB0;     % Number of AB individuals
            NBB = NBB0;     % Number of BB individuals
            
            % The simulation continues until only one genotype remains in the population

            while not((NAA ~= 0 && NAB == 0 && NBB == 0) || (NAA == 0 && NAB == 0 && NBB ~= 0))
                
                % Compute the probabilities of picking one AA, AB, or BB
                % individidual
                
                bAA = (1+delta)*NAA/((1+delta)*NAA+(1+h*delta)*NAB+NBB);
                bAB = (1+h*delta)*NAB/((1+delta)*NAA+(1+h*delta)*NAB+NBB);
                bBB = NBB/((1+delta)*NAA+(1+h*delta)*NAB+NBB);
                
                % Compute the probabilities of picking one A or B gamete
                
                bA = bAA+bAB/2;
                bB = bBB+bAB/2;
                
                % Compute the probabilities that the offspring is AA, AB,
                % or AB
                
                bAA_off = bA^2;
                bAB_off = 2*bA*bB;
                bBB_off = bB^2;

                % Generate a new offspring at random

                r1 = rand;
                
                if r1 <= bAA_off/(bAA_off+bAB_off+bBB_off)
                    
                    % Pick randomly one individual for death
                
                    r2 = rand;
                
                    if r2 <= NAA/(NAA+NAB+NBB)
                    
                        NAA = NAA-1;
                    
                    elseif r2 > NAA/(NAA+NAB+NBB) && r2 <= (NAA+NAB)/(NAA+NAB+NBB)
                    
                        NAB = NAB-1;
                    
                    elseif r2 > (NAA+NAB)/(NAA+NAB+NBB) && r2 <= 1
                        
                        NBB = NBB-1;
                    
                    end
                    
                    % A new offspring is born
                    
                    NAA = NAA+1;
                    
                elseif r1 > bAA_off/(bAA_off+bAB_off+bBB_off) && r1 <= (bAA_off+bAB_off)/(bAA_off+bAB_off+bBB_off)
                    
                    % Pick randomly one individual for death
                
                    r2 = rand;
                
                    if r2 <= NAA/(NAA+NAB+NBB)
                    
                        NAA = NAA-1;
                    
                    elseif r2 > NAA/(NAA+NAB+NBB) && r2 <= (NAA+NAB)/(NAA+NAB+NBB)
                    
                        NAB = NAB-1;
                    
                    elseif r2 > (NAA+NAB)/(NAA+NAB+NBB) && r2 <= 1
                    
                        NBB = NBB-1;
                    
                    end
                    
                    % A new offspring is born
                    
                    NAB = NAB+1;
                    
                elseif r1 > (bAA_off+bAB_off)/(bAA_off+bAB_off+bBB_off) && r1 <= 1
                    
                    % Pick randomly one individual for death
                
                    r2 = rand;
                
                    if r2 <= NAA/(NAA+NAB+NBB)
                    
                        NAA = NAA-1;
                    
                    elseif r2 > NAA/(NAA+NAB+NBB) && r2 <= (NAA+NAB)/(NAA+NAB+NBB)
                    
                        NAB = NAB-1;
                    
                    elseif r2 > (NAA+NAB)/(NAA+NAB+NBB) && r2 <= 1
                    
                        NBB = NBB-1;
                    
                    end
                    
                    % A new offspring is born
                    
                    NBB = NBB+1;
                    
                end

            end
            
            % Check whether the mutation took over the population
            
            if NAA ~= 0
                
                pfix = pfix+1;
                
            end
        
        end
        
        % Calculate the fixation probability
        
        pfix = pfix/Nrep;
   
end
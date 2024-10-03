clear;
global fitness trans_case birth_case death_case N group_size_male group_size_female number_state_male number_state_female number_state x_state y_state data;
get_fixation_probability(0.01, -2, -2, 5, 5)

function[fixation_probability] = get_fixation_probability(delta, h1, h2, N1, N2)
    % Predefined state transition vectors
    trans_case = [-1, 0, 0, 0; 
                  -1, 1, 0, 0;
                  0, -1, 0, 0;
                  0, 1, 0, 0;
                  1, -1, 0, 0;
                  1, 0, 0, 0;
                  0, 0, -1, 0;
                  0, 0, -1, 1;
                  0, 0, 0, -1;
                  0, 0, 0, 1;
                  0, 0, 1, -1;
                  0, 0, 1, 0;
                 ];
    birth_case = [2, 1, 2, 1, 0, 0];
    death_case = [0, 0, 1, 2, 1, 2];

    % Fitness matrix
    fitness = [1.0 + delta, 1.0 + h1 * delta, 1.0; 1.0 + delta, 1.0 + h2 * delta, 1.0];

    % Number of population states
    group_size_male =  N1;
    group_size_female = N2;
    number_state_male = (group_size_male + 2) * (group_size_male + 1) / 2;
    number_state_female = (group_size_female + 2) * (group_size_female + 1) / 2;
    number_state = number_state_male * number_state_female;
    x_state = ones(1, 12 * number_state);
    y_state = ones(1, 12 * number_state);
    data = zeros(1, 12 * number_state);

    % Solving fixation probability according to Markov transition matrix
    A = get_markov_matrix(group_size_male, group_size_female);
    A_new = A(2:number_state-1, 2:number_state-1);
    b = -A(2:number_state-1, number_state);
    x = bicgstabl(A_new, b, 1e-8, 1000);
    fixation_probability = x(1);
end

% Birth probability of three genotypes
function[birth_probability] = get_birth_probability(x1, x2, y1, y2, N1, N2)
    global fitness;
    weight_fitness_1 = times(fitness(1, :), [x1, x2, N1 - x1 - x2]);
    birth_probability_1 = [weight_fitness_1(1) + weight_fitness_1(2) / 2, weight_fitness_1(3) + weight_fitness_1(2) / 2];
    birth_probability_1 = birth_probability_1 / sum(birth_probability_1);
    weight_fitness_2 = times(fitness(2, :), [y1, y2, N2 - y1 - y2]);
    birth_probability_2 = [weight_fitness_2(1) + weight_fitness_2(2) / 2, weight_fitness_2(3) + weight_fitness_2(2) / 2];
    birth_probability_2 = birth_probability_2 / sum(birth_probability_2);
    birth_probability_AA = birth_probability_1(1) * birth_probability_2(1);
    birth_probability_BB = birth_probability_1(2) * birth_probability_2(2);
    birth_probability_AB = 1 - birth_probability_AA - birth_probability_BB;
    birth_probability = [birth_probability_AA, birth_probability_AB, birth_probability_BB] / 2;
end

% Death probability of three genotypes
function[death_probability] = get_death_probability(x1, x2, y1, y2, N1, N2, sex)
    if sex == 0
        death_probability = [x1, x2, N1 - x1 - x2];
    else
        death_probability = [y1, y2, N2 - y1 - y2];
    end
    death_probability = death_probability / sum(death_probability);
end

% Transition probability
function[transition_probability] = get_transition_probability(x1, x2, y1, y2, N1, N2, k)
    global birth_case death_case;
    sex = floor((k - 1) / 6);
    get_type = mod((k - 1), 6);
    birth_type = birth_case(get_type + 1);
    death_type = death_case(get_type + 1);
    birth_probability = get_birth_probability(x1, x2, y1, y2, N1, N2);
    death_probability = get_death_probability(x1, x2, y1, y2, N1, N2, sex);
    birth_probability = birth_probability(birth_type + 1);
    death_probability = death_probability(death_type + 1);
    transition_probability = birth_probability * death_probability;
end

% State mapping to number
function[index] = get_group_index(x1, x2, N)
    index = 0;
    space = N + 1;
    for t = 0:x1 - 1
        index = index + space;
        space = space - 1;
    end
    index = index + x2;
end

function[index] = get_index(x1, x2, y1, y2, N1, N2)
    global number_state_female;
    index_1 = get_group_index(x1, x2, N1);
    index_2 = get_group_index(y1, y2, N2);
    index = index_1 * number_state_female + index_2;
end

% Traverse the state pairs and construct the Markov transition matrix
function[markov_matrix] = get_markov_matrix(N1, N2)
    t = 1;
    global trans_case group_size_male group_size_female  number_state x_state y_state data;
    for i1 = 0:N1
        for i2 = 0:N1 - i1
            for j1 = 0:N2
                for j2 = 0:N2-j1                   
                    index_1 = get_index(i1, i2, j1, j2, group_size_male, group_size_female);
                    if index_1 == 0 || index_1 == number_state - 1
                        continue
                    end
                    change_probability = 0;
                    for k = 1:12
                        current_state = [i1, i2, j1, j2];
                        target_state = current_state + trans_case(k, :);
                        if target_state(1) < 0 || target_state(2) < 0 || target_state(3) < 0 || target_state(4) < 0 || target_state(1) + target_state(2) > N1 || target_state(3) + target_state(4) > N2
                            continue
                        end
                        index_2 = get_index(target_state(1), target_state(2), target_state(3), target_state(4), group_size_male, group_size_female);
                        transition_probability = get_transition_probability(i1, i2, j1, j2, group_size_male, group_size_female, k);
                        x_state(t) = index_1 + 1;
                        y_state(t) = index_2 + 1;
                        data(t) = transition_probability;
                        change_probability = change_probability + transition_probability;
                        t = t + 1;
                    end
                    x_state(t) = index_1 + 1;
                    y_state(t) = index_1 + 1;
                    data(t) = -change_probability;
                    t = t + 1;
                end
            end
        end
    end
    markov_matrix = sparse(x_state, y_state, data, number_state, number_state);
end

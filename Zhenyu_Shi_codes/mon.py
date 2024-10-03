import time
import numpy
import numpy as np
import datetime


# Predefined state transition vectors
case = np.array([[-1, 0],
                [-1, 1],
                [0, -1],
                [0, 1],
                [1, -1],
                [1, 0],
                ]
                )
birth_case = np.array([2, 1, 2, 1, 0, 0])
death_case = np.array([0, 0, 1, 2, 1, 2])
# Number of population states
group_size = 10
number_state = (group_size + 2) * (group_size + 1) // 2
markov_matrix = np.zeros((number_state, number_state), dtype=np.float32)


# Birth probability of three genotypes
def get_birth_probability(x1, x2, N):
    weight_fitness_1 = fitness * np.array([x1, x2, N - x1 - x2])
    birth_probability_1 = np.array([weight_fitness_1[0] + weight_fitness_1[1] / 2, weight_fitness_1[2] + weight_fitness_1[1] / 2])
    birth_probability_1 /= birth_probability_1.sum()
    weight_fitness_2 = fitness * np.array([x1, x2, N - x1 - x2])
    birth_probability_2 = np.array([weight_fitness_2[0] + weight_fitness_2[1] / 2, weight_fitness_2[2] + weight_fitness_2[1] / 2])
    birth_probability_2 /= birth_probability_2.sum()
    birth_probability_AA = birth_probability_1[0] * birth_probability_2[0]
    birth_probability_BB = birth_probability_1[1] * birth_probability_2[1]
    birth_probability_AB = 1 - birth_probability_AA - birth_probability_BB
    birth_probability = np.array([birth_probability_AA, birth_probability_AB, birth_probability_BB])
    return birth_probability / 2


# Death probability of three genotypes
def get_death_probability(x1, x2, N):
    death_probability = np.array([x1, x2, N - x1 - x2], dtype=float)
    death_probability /= death_probability.sum()
    return death_probability


# Transition probability
def get_transition_probability(x1, x2, N, k):
    birth_type = birth_case[k]
    death_type = death_case[k]
    birth_probability = get_birth_probability(x1, x2, N)[birth_type]
    death_probability = get_death_probability(x1, x2, N)[death_type]
    return birth_probability * death_probability


# State mapping to number
def get_index(x1, x2, N):
    index = 0
    space = N + 1
    for i in range(x1):
        index += space
        space -= 1
    index += x2
    return index


# Traverse the state pairs and construct the Markov transition matrix
def get_markov_matrix(N):
    for i1 in range(N + 1):
        for i2 in range(N + 1 - i1):
            index_1 = get_index(i1, i2, group_size)
            change_probability = 0
            for k in range(6):
                current_state = np.array([i1, i2])
                target_state = current_state + case[k]
                if target_state.max() > N or target_state.min() < 0 or target_state[0] + target_state[1] > N:
                    continue
                index_2 = get_index(target_state[0], target_state[1], group_size)
                transition_probability = get_transition_probability(i1, i2, group_size, k)
                markov_matrix[index_1, index_2] = transition_probability
                change_probability += transition_probability
            markov_matrix[index_1, index_1] = 1 - change_probability


fitness = np.array(np.array([1 + 0.01, 1 - 0.01, 1.0]))
group_size = 10
number_state = (group_size + 2) * (group_size + 1) // 2
get_markov_matrix(group_size)
# Solving fixation probability according to Markov transition matrix
markov_matrix_new = markov_matrix[1:number_state-1, 1:number_state-1]
b = markov_matrix[1:number_state-1, -1]
a = np.linalg.solve(markov_matrix_new - np.eye(number_state - 2), -b)


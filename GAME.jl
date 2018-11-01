

A1 = "AAGTAAGAAAATGAAGCCCCTTGTCCAGCAGTCAGCTTAAGGCGTGGTGTAGATAATATATGCATTAGACCTTAGTTAAGTCACGCAGTCACAATTAGGCCGAATTAcCGGGCCAGGTTCCCACGTTAGGGACCCTGCATTGCTCGTTTAGGATATGTCTTTTGCTGCGTTTGCGAGACACAAGGTAGCAAGGATACGTCCGGGCACTCTCTTACA"
A2 = "TACGTAGCTACTTGCAGGACCGCTTAAGTCTGCTTTGATGTCCTGGTTCGTGGCAGGCGGGGAGTCGAGATATCTCTAATTCGCCCGAACGGTTATAATTAGaCCGAATTAATACCAAGGCAGTCTGCCGGTGTTGTATGGTGTCTTGACGAGAATAGTAAGAGAATGAGGTGTCCCGGCGTGATTGGCTGGTGGGTGTTTTACTAAGAGCCTGAA"
A3 = "TCGCAAGCCTCCCCACTCCATCCGTTCAAGGTTCACCGACGTCCTAACTGCCGTTCCCAATTAGGCgGAATTAACTGGCTGGGGAGTGGGATTGAAACCAGGCTTTTTATCAAGAGGTACGGCTGTCTGTAGCGTGACGTGAGGCTCAGCGGCTTAACATAGACACCGTCTGCACATAATTTCTCGCTTCATGGTCTCGAAATTGACAAATGATCC"
A4 = "GCACAGCTGGCTCACCTTGTATGCACAGGCTGACTTCATGGCACGCCGATTCCAGCCGGGAATaAGGaCGAATTAACGTCTACAAATAGCTGAGAGACCTTCTAACATAAGAAGCGCTGTAGTCGCTGCCCCAACGTAGATTGGCGGGTCCAACACTAGGCATGTTGAACCCAATTTACGGCGGCTCGCCCGTTACAATCGTAGAGCCATACCGCC"

### population = Array(Int64,12,2)
populationSize = 10
seq = [A1,A2,A3,A4]
population = hcat([rand(1:length(seq[1]),length(seq)) for i in 1:populationSize]...)
population = transpose(population)
### seq = hcat([seq[i] for i in 1:length(seq)])


############################## Crossover function code ###########################################################################################
import Random
Random.seed!(1234)

function crossover(initial_population)    
    shuffled_initial_population = initial_population[Random.randperm(length(initial_population))] # Shuffles/randomizes initial population
    crossover_population = [] # Will hold children produced by crossover (N/2 pairs)
    
    for i in 1:2:length(shuffled_initial_population) # index by 2
        if(i+1 <= length(shuffled_initial_population))
            A = shuffled_initial_population[i]
            B = shuffled_initial_population[i+1]
            crossOverPoint = rand(1:length(A)) # cross over point; we assume A and B are of the same length (m)
        
            A_prime = []
            B_prime = []
            for j in 1:crossOverPoint
                push!(A_prime, A[j])
                push!(B_prime, B[j])
            end
            for j in crossOverPoint+1:length(A)
                push!(A_prime, B[j])
                push!(B_prime, A[j])
            end

            push!(crossover_population, A_prime)
            push!(crossover_population, B_prime)
        end
    end
    
    return crossover_population
end
  
##############################################################################################################################################


###### Helper functions to generate sample data ##############################################################################################
import Random

# generate random set of DNA sequences of given size and quantity
function generate_random_sequences(number_of_sequences, sequence_length)
    generated_sequences = []
    for i in 1:number_of_sequences
        generated_sequence = Random.randstring("ACGT", sequence_length)
        push!(generated_sequences, generated_sequence)
    end
    return generated_sequences
end

# generate a random population of potential motif locations of a given size for a given set of DNA sequences
function generate_random_initial_population(population_size, number_of_sequences, sequence_length, motif_width)
    generated_population = []
    # seeds = []
    for i in 1:population_size
        motif_start_point_range = 1:sequence_length - motif_width
        member_size = number_of_sequences
    	  # seed = rand(blah)
    	  # Random.seed!(seed)
    	  # push!(seeds, seed)
        generated_member = rand(motif_start_point_range, member_size) # Generate a member of the population of size t
        push!(generated_population, generated_member)
    end
    return generated_population
end

Random.seed!(1234);

generate_random_initial_population(3,5,10,4)


##############################################################################################################################################

###################### Mutation functions ####################################################################################################

# IP: mutation_probability:Float(0,1)
# OP: boolean value indicating whether to mutate a single bit of a candiate in the population
function to_mutate(mutation_probability)
    number = rand(1)[1]
    if number  <= mutation_probability
        return true
    else
        return false
    end
end

# Test function
# function test_mutation_occurences():
#     P = 0.000000000001
#     for i in 1:50000000
#         if to_mutate(P)
#             println(true)
#     end
# end
  
# IP: initial_population: List<List:Int>, mutation_probability:Float(0,1), sequence_length:Int
# OP: mutated_population: List<List:Int>
# Algorithm: Perform mutations over a population of candidate motifs with a given mutation probability assuming uniform sequence legnths
function mutation(initial_population, mutation_probability, sequence_length)    
    mutated_population = [] # mutated_population
    for i in 1:length(initial_population)
        current_candidate = deepcopy(initial_population[i])
        for j in 1:length(current_candidate)
            if to_mutate(mutation_probability)
                current_candidate[j] = rand(1:sequence_length)
            end
        end
        push!(mutated_population, current_candidate)
    end
    return mutated_population
end


# Test function
function test_mutations()
    sequence_length = 200
    population_size = 10
    number_of_sequences = 200
    initial_population = generate_random_initial_population(population_size, number_of_sequences, sequence_length)
    mutated_population = mutation(initial_population, 0.001, sequence_length)

    # println(typeof(initial_population), typeof(mutated_population))

    for (i, j) in zip(initial_population, mutated_population)
        if any(i .!= j)
            println(i[i .!= j])
            println(j[i .!= j])
            println("*********************************")
        end
    end
    return initial_population, mutated_population
end

# initial_population, mutated_population = test_mutations()

### Helper functions to analyze mutation results ###

# Count how many members of the population were mutated
function count_mutated_members(initial_population, mutated_population)
    mutated_members = 0
    for (i, j) in zip(initial_population, mutated_population) 
        if any(i .!= j)
            mutated_members += 1
        end
    end
    return mutated_members
end

import StatsBase

# Count how many point mutations occurred across all members of the population
function count_mutations(initial_population, mutated_population)
    mutations = 0
    for (i, j) in zip(initial_population, mutated_population)
        mask = (i .!= j)
        mutations += StatsBase.countmap(mask)[true]
    end
    return mutations
end



############################## tournament Selection function code ###########################################################################################

function Score(Alignment)
    k = length(Alignment[1])
    profile = Dict('A' => [0 for i in 1:k],
            'C' => [0 for i in 1:k],
            'T' => [0 for i in 1:k],
            'G' => [0 for i in 1:k]) ### 1 profile
    
    for i in 1:length(Alignment) ### For each kmer
        for j in 1:length(Alignment[i]) ### traverse each base
            base = uppercase(Alignment[i][j])
            profile[base][j] += 1
        end
    end
    m = sum(values(profile))[1]
    score = 0
    for col in 1:length(profile['A'])
        score += (m - max(profile['A'][col], profile['C'][col], profile['T'][col], profile['G'][col]))
    end
    return(score)
end

import Random

function profile(population,seq,width)
    p = [Dict('A' => 0, 'C' => 0, 'T' => 0, 'G' => 0) for i in 1:width] ### build a profile
    parray = [p for i in 1:size(population)[1]] ### build an array of profiles
    for i in 1:size(population)[1] ### For each individual
        for j in 1:size(population)[2] ### traverse each seq
            if (population[i,j]+width > length(seq[j]))
                br = length(seq[j])
            else
                br = population[i,j]+width
            end
            for c in 1:length(seq[j][population[i,j]:br-1]) ### Traverse each letter
                parray[i][c][uppercase(seq[j][c+population[i,j]])] += 1
            end
        end
    end
    return(parray)
end

function tournamentSelection(population,tSize,seq,motifWidth)
    ### Population is group into tournaments of size tSize.
    ### Output Winners of tournaments
    
    tournamentPopulation = population[Random.randperm(size(population)[1]),:]
    Winners = []
    
    for j in 1:tSize:(size(tournamentPopulation)[1]-tSize+1)
        println(j)
        println(tournamentPopulation[j:j+tSize-1,:])
        bracket = tournamentPopulation[j:j+tSize-1,:]
        
        topScore = typemax(Int64) ### Worst possible score is large number
        bracketWinner = -1
        for b in 1:size(bracket)[1]
            CurrentAlignment = [seq[idx][bracket[b][idx]:bracket[b][idx]+motifWidth-1] for idx in 1:length(seq)]
            # println(CurrentAlignment)
            # println(Score(CurrentAlignment))
            if(topScore > Score(CurrentAlignment) )
                topScore = Score(CurrentAlignment)
                println("Found new top score:", topScore)
                bracketWinner = bracket[b]
            end
        end
       push!(Winners,bracketWinner)
    end
    return Winners
end
  
  












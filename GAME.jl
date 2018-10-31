



function tournamentSelection(x::AbstractString)
    
end

tournamentPopulation = hcat([population[:,Random.randperm(size(population)[2])] for i in 1:tSize]...) ### Self Tournaments are possible here.

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

############################## tournament Selection function code ###########################################################################################
import Random
q1 = 0
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

function tournamentSelection(population,tSize,seq)
    ### Population is group into tournaments of size tSize.
    ### Output Winners of tournaments
    
    tournamentPopulation = population[Random.randperm(size(population)[1]),:]

    for j in 1:tSize:size(tournamentPopulation)[2]-tSize
        println(j)
        p = profile(tournamentPopulation[j:(j+tSize-1),:],seq,width)
    end
    return crossover_population
end
  


  
##############################################################################################################################################



### HI


x = 5
print(x)


# ADJUST:  A(a1,…,ai,…,am)  
# i ? 1 
# DO 
#      a'i=argmax  ?ent(a1,?…?,?ai,?…?,?am), 1 = ai = li - w + 1 
#      ai?a'i  
#     i ? i + 1 
#     if i > m then i ? 1 
# UNTIL no further improvements obtained 
  




# SHIFT:  (A(a1,…,ai,…,am))  
# k' = argmax  ?ent(a1+k,………,am+k),-w=k=w  
# A?A(a1+k'………,am+k')  
# Note 1: if ai = 0 then ai + k = 0 for all k (no added sites) 
# Note 2: If ai + k < 0 or ai + k > li - w + 1, then set ai + k = 0 




BEGIN 
    Initialization: i ? 0 
        Setting parameters: 
            population size N = 500 
            mutation rate r = 0.001 
            maximum generation G = 3000: 
        Generating initial population P0: 
    Repeat: i ? i + 1 
        Mutate individuals;  
        Crossover individuals;  
        Selection of individuals; 
    Until (i = G or convergence) 
    Choose best individual Aopt; 
    Repeat 
        ADJUST(Aopt); 
        SHIFT(Aopt); 
    Until no further improvements obtained 
    PWM-scan on Aopt to extract additional weaker motif sites 
END 



  
A1 = AAGTAAGAAAATGAAGCCCCTTGTCCAGCAGTCAGCTTAAGGCGTGGTGTAGATAATATATGCATTAGACCTTAGTTAAGTCACGCAGTCACAATTAGGCCGAATTAcCGGGCCAGGTTCCCACGTTAGGGACCCTGCATTGCTCGTTTAGGATATGTCTTTTGCTGCGTTTGCGAGACACAAGGTAGCAAGGATACGTCCGGGCACTCTCTTACA
A2 = TACGTAGCTACTTGCAGGACCGCTTAAGTCTGCTTTGATGTCCTGGTTCGTGGCAGGCGGGGAGTCGAGATATCTCTAATTCGCCCGAACGGTTATAATTAGaCCGAATTAATACCAAGGCAGTCTGCCGGTGTTGTATGGTGTCTTGACGAGAATAGTAAGAGAATGAGGTGTCCCGGCGTGATTGGCTGGTGGGTGTTTTACTAAGAGCCTGAA
A3 = TCGCAAGCCTCCCCACTCCATCCGTTCAAGGTTCACCGACGTCCTAACTGCCGTTCCCAATTAGGCgGAATTAACTGGCTGGGGAGTGGGATTGAAACCAGGCTTTTTATCAAGAGGTACGGCTGTCTGTAGCGTGACGTGAGGCTCAGCGGCTTAACATAGACACCGTCTGCACATAATTTCTCGCTTCATGGTCTCGAAATTGACAAATGATCC
A4 = GCACAGCTGGCTCACCTTGTATGCACAGGCTGACTTCATGGCACGCCGATTCCAGCCGGGAATaAGGaCGAATTAACGTCTACAAATAGCTGAGAGACCTTCTAACATAAGAAGCGCTGTAGTCGCTGCCCCAACGTAGATTGGCGGGTCCAACACTAGGCATGTTGAACCCAATTTACGGCGGCTCGCCCGTTACAATCGTAGAGCCATACCGCC

strand_length = 216
motif_length = 5
  

M1 = [5, 6, 7, 8]
M2 = [2, 4, 10, 14]
M3 = [3, 9, 27, 28]
M4 = [35, 46, 47, 58]
M5 = [33, 42, 41, 23]
M6 = [44, 50, 51, 32]

list1 = [M1, M2, M3, M4, M5, M6]
list2 = [M2, M3, M5, M1, M6, M4]

1:length(initial_population)
derrangement(1:initial_population)


  
initial_population = [M1, M2, M3, M4, M5]

initial_population = [(5), (5), (5), (5)] ### A = Array(Float64,1,1) ### Init as matrix 

https://rosettacode.org/wiki/Permutations/Derangements#Julia

after_crossover_population = append!(initial_population, initial_population)


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

  julia> Random.seed!(1234);

julia> x2 = rand(2)
2-element Array{Float64,1}:
 0.590845
 0.766797
  
  generate_random_initial_population(3,5,10,4)
3-element Array{Any,1}:
 [2, 6, 1, 2, 4]
 [6, 4, 2, 6, 4]
 [5, 2, 3, 1, 3]

#Mayank Op:
[5, 4, 4, 5, 4]
[3, 2, 4, 3, 5]
[6, 3, 4, 1, 3]

  
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


##############################################################################################################################################
  ### Ian's sad work

  ### This is not done yet. I struggled for a while with vectorization, like this: a = [1,2,3]; sin.(a); It seems that you guys are using arrays of arrays or arrays of strings
  ### rather than matrices. It will be fine either way, but we should probably make sure everyone is on the same page with that. Array{Float64,1}
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

function tournamentSelection(population,tSize,seq)
    ### Population is group into tournaments of size tSize.
    ### Output Winners of tournaments
    
    tournamentPopulation = population[Random.randperm(size(population)[1]),:]

    for j in 1:tSize:size(tournamentPopulation)[2]-tSize
        println(j)
        p = profile(tournamentPopulation[j:(j+tSize-1),:],seq,width)
    end
    return crossover_population
end
  
function score()
  ### I was supposed to write this and I didn't because data structures are hard! #sad
  
end











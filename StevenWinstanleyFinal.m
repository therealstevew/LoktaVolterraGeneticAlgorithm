%Steven Winstanley Sheridan College S/N: 991457192
%Clear Screen
clear;
%Loop count
k = 1;
%Generate Intial Population of x amount of species 
population = genisis(2);
%Storage for first
pop1 = {};
%Storage for second
pop2 = {};
average1 = 0;
average2 = 0;
%Main Loop
while k < 100
    %Check if only one species remains
    if(length(population) == 1)
        %check if it belongs to population 1
        if length(population{1}) == length(pop1)
            %Set last value for plotting
            x1(k) = 0;
            t1(k) = k;
            t(k) = k;
            x(k) = length(pop1);
            disp("Population 1 survies");
        %check if belongs to pop 2
        else
            x(k) = 0;
            t(k) = k;
            t1(k) = k;
            x1(k) = length(pop2);
            disp("Population 2 survies");
        end
        break;
    end
    %Store this instance of each population
    pop1 = population{1};
    pop2 = population{2};
    disp(pop1);
    disp(pop2);
    %Save population data for plotting
    x(k) = length(pop1);
    x1(k) = length(pop2);
    t(k) = k;
    t1(k) = k;
    if(length(pop1) ~= 0)
        average1 = average1 + sum([pop1{:}]);
    end
    if length(pop2) ~= 0 
        average2 = average2 + sum([pop2{:}]);
    end
    %advance counter
    k = k + 1;
    %Create next generation 
    population = make_next_gen(population);
end
%plot values
disp("Average Fitness for Population 1");
disp(average1/k);
disp("Average Fitness for Population 2");
disp(average2/k);

plot(t, x);
hold on;
plot(t1, x1);
xlabel("Time (ecological cycles)");
ylabel("Population of Species");
legend('Population1', 'Population 2')
%Inital Creation function
function out = genisis(num_of_species)
    %Declare population
    universe = {};
    %Create x amount of species
    for i = 1:num_of_species
        species = {};
        %Generate individuals for each species random between 5 - 15
        for j = 1:randperm(15, 5)
            %Generate individual value between 1-4
            species{end+1} = randperm(4, 1);
        end
        %add newly created species to overall population
        universe{end+1} = species;
    end
    %return population
    out = universe;
end
%Calculate pressure
function out = pressure(direct, indirect)
    %Get species in question pressure on the enviroment 
    ax = sum([direct{:}]);
    %init ay
    ay = 0;
    pop_size = 0;
    %get other species and sum their total effect
    for i = 1:length(indirect)
        other_pop = indirect{i};
        for j = 1:length(other_pop)
            ay = ay + other_pop{j};
        end
        pop_size = pop_size + length(indirect{i});
    end
    ax = ax/ ay;
    ax = ax * pop_size;
    out = ax;
end
%Calculating growth rate
function out = calc_growth_rate(size, pressure, time)
    %Enviromental Carrying capacity 
    k = 20;
    %rate of growth
    r = 5;
    %init growth Rate
    growth_rate = -1;
    %check the size of population
    if(size ~= 0)
        growth_rate = ((r*size)*((k - size) - pressure) / size) * time;
    end
    %return growth rate
    out = growth_rate;
end
%Calculating the fitness of a population
function out = calc_fitness(population)
    % sum values of population
    out = sum([population{:}])
end
%This is used to select parents for offspring
function out = choice_by_roulete(population, fitness_sum)
    %init offset and normal for roulette calculation
    offset = 0;
    normal = fitness_sum;
    %convert given cell to mat and sort
    lowest_fitness = cell2mat(population);
    lowest_fitness = sort(lowest_fitness);
    %assign the lowest fitness
    lowest_fitness = lowest_fitness(end);
    size = length(population);
    %check if there is a need to offset the current pop within the wheel
    if(lowest_fitness < 5)
        offset = -lowest_fitness;
        normal = normal + offset * size;
    end
    
    %generate choice
    draw = rand();
    accum = 0;
    %select pop from population
    for i = 1:length(population)
        fitness = population{i} + offset;
        prob = fitness / normal;
        accum = accum + prob;
        if (draw < accum)
            out = population{i};
            break;
        end
        
    end    
end
%Creating the next generation of pops given growth rate
function out = make_next_gen(population)
    %Loop through every species in the population
    for i = 1:length(population)
        fitness_sum = 0;
        try
            species = population{i};
        catch
            
        end
        %Break if one species remains 
        if (length(population) == 1)
            break;
        end
        %Make a species go extinct if there is only 1 member left
        if (length(species) == 1)
            population(i) = [];
            continue;
        end
        %create a population without species i to calculate the presssure
        %from other species
        other_pop = population;
        other_pop{i} = [];
        %Calculate the growth rate of species
        pop_growth = calc_growth_rate(length(species), pressure(species, other_pop), 0.175);
        pop_growth = round(pop_growth);
        sorted_by_fitness = species;
        %calculate fitness of species i
        fitness_sum = calc_fitness(species);
        % there cant be stagnat population so generating a random number to
        % simulate variance
        if(pop_growth == 0)
            pop_growth = pop_growth + randperm(5, 1);
        end
        %if a population is decreasing in number
        if(pop_growth < 0)
           %Calculate the effect of growth rate on current species 
            dif = length(species) - pop_growth * -1;
            %if the difference between current pop and growth is > 0 kill
            %off the species
            if (dif <= 0)
                population(i) = [];
                continue;
            else
                %if the species survives after negative growth reduce that
                %many from the current population
                for j = 1:(pop_growth * -1)
                    species(end) = [];
                end
            end
        else
            % if pop growth is positive create new pops 
            for j = 1:pop_growth
                %Select parentsd 
                first_choice = choice_by_roulete(sorted_by_fitness, fitness_sum);
                second_choice =  choice_by_roulete(sorted_by_fitness, fitness_sum);
                %give new offspring parents genes
                individual = crossover(first_choice, second_choice);
                %calculate if there is a mutation within the child 
                if (rand() < 0.5)
                    individual = mut(individual);
                end
                %Add new offspring to current species
                species{end+1} = individual;
            end
        end
        %Append altered species to population
        population{i} = species;
    end
    out = population;
end

%Offspring function
function out = crossover(a, b)
    out = (a/2) + (b/2); 
end
%Mutation function
function out = mut(individual)
    child = individual;
    child = child + randperm(2, 1);
    out = child;
end

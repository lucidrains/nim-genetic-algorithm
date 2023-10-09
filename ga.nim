assert NimMajor >= 2, "must be nim 2.0"

# imports

import std/[math, random, terminal, strformat, algorithm, sequtils, sugar]

# init

randomize()

# gene

type
  Gene = object
    code: string
    cost: int = 1e7.int
  GeneRef = ref Gene

proc randomize(gene: GeneRef, length: int) =
  assert length > 2

  while len(gene.code) < length:
    gene.code = gene.code & chr(rand(255))

# finds a random index in the code and adds or subtracts 1 from the ord randomly

proc mutate(gene: GeneRef) =
  let rand_index = rand(gene.code.len - 1)
  let mutation = (if rand(1.0) < 0.5: 1 else: -1)

  var code_seq_char = cast[seq[char]](gene.code)
  let char_at_index = code_seq_char[rand_index]

  var char_code = ord(char_at_index) + mutation
  char_code = min(max(char_code, 0), 255)

  code_seq_char[rand_index] = chr(char_code)
  gene.code = cast[string](code_seq_char)

# recombine genetic codes

proc mate(gene1, gene2: GeneRef): GeneRef = 
  assert gene1.code.len == gene2.code.len
  let pivot = rand(1..(gene1.code.len - 2))
  let new_code = gene1.code[0..<pivot] & gene2.code[pivot..^1]
  return GeneRef(code: new_code)

proc calc_cost(gene: GeneRef, goal: string) =
  assert len(gene.code) == len(goal)

  var cost: int = 0

  for i in 0..<len(goal):
    let diff: int = ord(goal[i]) - ord(gene.code[i])
    cost += pow(diff.float, 2).int

  gene.cost = cost

proc cmp_genes(gene1, gene2: GeneRef): int =
  gene1.cost - gene2.cost

# population

type
  Population = object
    size: int
    pool: seq[GeneRef]
    goal: string 
    generation: int
    keep_fittest_frac: float = 1.0
    mutate_prob: float = 0.5
    solved: bool = false
  PopulationRef = ref Population

proc init_population(pop: PopulationRef) =
  for _ in 0..<pop.size:
    let gene = GeneRef()
    randomize(gene, len(pop.goal))
    pop.pool.add(gene)

proc sort_by_fitness(pop: PopulationRef) =
  pop.pool.sort(cmp_genes, Ascending)

proc keep_fittest(pop: PopulationRef) =
  pop.sort_by_fitness()
  let num_survive: int = (pop.size.float * pop.keep_fittest_frac).float.int
  pop.pool = pop.pool[0..<num_survive]

proc breed_fittest(pop: PopulationRef) =
  assert pop.pool.len >= 2

  while len(pop.pool) < pop.size:
    pop.pool.shuffle()
    let child = mate(pop.pool[0], pop.pool[1])
    assert child.code.len == pop.goal.len
    pop.pool.add(child)

# next generation

proc next_generation(pop: PopulationRef) =

  # if population not started, initialize with random genes

  if pop.generation == 0:
    init_population(pop)

  # calculate fitness

  for gene in pop.pool:
    gene.calc_cost(pop.goal)

  # keep only the fittest

  pop.keep_fittest()

  # breed the fittest

  pop.breed_fittest()

  # mutate

  for gene in pop.pool:
    if rand(1.0) < pop.mutate_prob:
      gene.mutate()

    gene.calc_cost(pop.goal)

  pop.solved = pop.pool.any(g => g.cost == 0)

  if pop.solved:
    pop.sort_by_fitness()

  pop.generation += 1

# display population of genes

proc display(pop: PopulationRef) =
  erase_screen()

  echo fmt"generation: {pop.generation}", "\n"

  for gene in pop.pool:
    echo fmt"{gene.code} ({gene.cost})"

  echo "\n"

proc main() = 
  # create a gene population

  let population = PopulationRef(
    goal: "Attention is all you need",
    size: 20,
    keep_fittest_frac: 0.25,
    mutate_prob: 0.5
  )

  # while not solved, do another generation

  while not population.solved:
    population.next_generation()
    population.display()

when is_main_module:
  main()

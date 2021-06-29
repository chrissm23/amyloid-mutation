# amyloid-mutation
Setup to run free energy calculations on amyloid mutations

I want to provide a workflow to automatically setup all or at least most things required for free energy calculations on amyloid mutations.
The project contains python and shell files to achieve this objective which can hopefully be merged afterwards into one single program.

The steps of the workflow are:

1. **Get data about the structure**
  * What residues are present in the stucture?
  * How many chains are there?
  * Arrange how the chains will be mutated.
  * Create amyloid mutations.
  * What is the charge of WT and mutants?
  * Determine what the best order of adding ions is:
    * If they have same charge sign, add ions to least charged first.
    * If they have oppposite charge sign, add ions to the one that has the same charge sign as WT first, then add ions to other mutant, then balance charge of second solvation to get -charge of first structure.

2. **Add water and ions**
  * Solvate WT and all mutants.
  * Strip water.
  * Add ions according to above order and strip water with ions and, just in case, amyloid with ions as well.

3. **Intercalate chains of each pair of structures and merge**
  * Create a PDB for each mutant intercalated with WT or previous mutant.
  * Get chains and residues to add to tiMerge.
  * Combine with respective solvent.
  * Merge according to above data.

4. **Copy and complete templates**
  * Copy tmeplates to corresponding directories and complete with corresponding masks.

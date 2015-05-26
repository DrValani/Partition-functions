Computing exact 2D/3D Partition Functions
===================

This repo contains my source code for computing exact partition functions on 2 and 3 dimensional lattices. 

Firstly: a warning on the code... this was written before I read Robert Martin's Book: Clean Code. Also before I gained  experience as a professional Software Developer. So expect a mess, there are no unit or integration tests. My testing was based on writing a brute force code to calculate results of smaller managable lattices sizes. Then with each iteration, the new (hopefully) faster program's results would be tested against the previous iteration. Plotting the zero's of each partition function on a 2d plane revealed errors too (more on that in my thesis [http://ethos.bl.uk/OrderDetails.do?uin=uk.bl.ethos.546769]). Not to mention a check against Pearson's exact result [http://journals.aps.org/prb/abstract/10.1103/PhysRevB.26.6285].


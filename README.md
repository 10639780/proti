# proti
Proteine folding tool

Als we spiegelbeelden met een constraint wegfilteren kunnen we de berekenbaarheid van het probleem verbeteren zonder een mogelijk optimale uitkomst weg te gooien. Werken i  

Wat betreft de harde constraint dat het molecuul zichzelf niet mag kruisen, kunnen we benaderen hoeveel mogelijkheden hiermee worden weggegooid? 

Okke heeft ons aangeraden om een verkenning te doen van heuristieke algorithmen, welke werkt daarbij het beste? Monte carlo is hiervan een eerste voorbeeld, maar het is goed om een indruk te krijgen van wat een andere aanpak met het probleem doet. Binnen het kader van dit vak hebben we niet de mogelijkheid om uitsluitend een exacte benadering te nemen.

mogelijke keuzes: 

ant colony 
iterative deepening
genetisch geinspireerd 
A*

Een mooi onderwerp voor de eindpresentatie is waar een heuristisch algorithme het van een exact algorithme moet overnemen. Welke variabelen zijn van invloed in deze afweging? 

Goed bijhouden van de dingen die we tot nu toe geleerd hebben aan de hand van monte carlo kan twee dingen opleveren. Ten eerste kunnen we aan de hand hiervan beter begrijpen of een ander algorithme beter presteert of niet. Ten tweede is de combinatie van elementen dan ook te beoordelen. Het kan interessant zijn om een heuristisch algorithme te bouwen dat andere elementen combineert, hiervoor moeten we de mogelijkheden verkennen. Een grafische weergave van de resultaten is hierbij een must. Code tonen wordt tijdens de eindpresentatie niet gewaardeerd. Pseudocode tonen daarentegen wel!

Okke heeft een functie uitlegt die op een exacte manier (door een grote overschatting te gebruiken) kan voorspellen of een bepaalde vertakking in de boom -door het potentiele aantal punten in die vertakking te nemen- een verbetering kan geven in de uiteindelijke score. Deze pruning methode gebeurt aan het einde van de boom structuur. Het verschil met de oppervlakte benadering is dat hij een stuk specifieker is en dus niet hoger in de boom toegepast kan worden. Toch is het een andere uitdrukking van dezelfde denkrichting en daarom interessant voor onze case. 

Er zijn nu dus twee richtingen die we gaan volgen. Ten eerste hebben we het exacte algorithme dat voorlopig prioriteit heeft. Het is op dit moment belangrijk om zoveel mogelijk mogelijkheden wat betreft pruning te vinden. Tegelijkertijd moeten we voor de volgende presentatie helder kunnen vertellen wat (annealing) monte carlo en andere heuristische algorithmen doen met de state space: waar zijn ze goed in? Waar zijn ze minder goed in? Op een gegeven moment zullen we van het exacte pad af moeten stappen en het is belangrijk dat we weten wat de mogelijkheden zijn als we eenmaal zover zijn. 

Het visualiseren van de voortgang zal tot en met de eindpresentatie van belang zijn. Kunnen we al data gaan verzamelen? Een presentatie is een verhaal van begin tot eind, wat kunnen we tot nu toe vertellen over het begin?
 
Overzicht: 
- Exact: pruning met de oppervlaktemethode
- Algemene tool voor visualisering, voor alle programma's hetzelfde -> uitcommenten wat je niet nodig hebt.
- Nieuwe heuristische benaderingen ant colony, iterative deepening, genetische methode.

termen

structuur = molecuul 

element = amino = H/P/C

HP model protein folding veel papers over!

atoom bestaat niet (we werken op een ander level)






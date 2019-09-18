# Endorse
Predikce vlastností EDZ s vlivem na bezpečnost a spolehlivost hlubinného úložiště radioaktivního odpadu.

##Cíl projektu 
Vytvořit SW nástroj a metodiku, pro predikci veličin charakterizujících bezpečnost dílčí části úložiště
(tzv. *indikátorů bezpečnosti*) na základě geofyzikálních měření. To zahrnuje:

1. Sestavení modelu transportu kontaminace skrze EDZ od (náhodných) úložných kontejnerů do hypotetické poruchy. 
Zahrnutí předpokládané geometrie úložiště s velikostí do 100m.
2. Definice vhodných indikátorů bezpečnosti jakožto veličin odvozených od výzledků modelu transportu.
3. Tvorbu menších modelů pro identifikaci parametrů transportního modelu na základě předpokládaných průzkumů 
a geofyzikálních měření.
4. Aplikaci vhodných stochastických výpočetních metod pro predikci rozdělení indikátorů bezpečnosti a parametrů 
transportního modelu se zahrnutím chyb měření a dalších podstatných neurčitostí použitých modelů

## Zdroje
[Přehled řešitelů](https://docs.google.com/document/d/1R8CBU9197brrruWGahVbE7_At2S2V51J6JV5bgs-kxQ/edit#heading=h.e1t1yg8nyvaz)
[repozitář Flow123d](https://github.com/flow123d/flow123d) 
- simulátor transportních a mechanických procesů v rozpukaném porézním prostředí
[repozitář MLMC](https://github.com/GeoMop/MLMC)
- metoda multilevel Monte Carlo, generování náhorných polí a puklinových sítí, 
  maximal entropy method pro rekonstrukci hustoty pravděpodobnosti
  

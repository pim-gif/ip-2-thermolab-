---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.19.6
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

# Labjournaal oploswarmte van zout
## Algemeen
In een jupyter labjournaal maak je gebruik van markdown. Voor veel gebruikte notaties, [check de documentatie](https://jupyter-notebook.readthedocs.io/en/stable/examples/Notebook/Working%20With%20Markdown%20Cells.html).


```python
#import necessary libraries
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit 
```

**Namen & studentnummer**

- Cem Aras 685667
- Pim Kroes 6806821   
- Jurre van Woerden 6805183


**Titel van het experiment:**  Oploswarmte van zout


**Datum:**  26/11/2026


**Doel van het experiment:**  
In dit experiment bepalen we de oploswarmte van zout in water.

Metingen van de verandering van de weerstand van een “sensor” worden vaak gebruikt voor de bepaling van veranderingen in fysische grootheden. Hierbij kunnen wij denken aan veranderingen in temperatuur, druk, uitrekking of verandering van de weerstand door een magneetveld. Weerstandsveranderingen kunnen met grote gevoeligheid gemeten worden. Met een 6½ digit DMM kan in principe een relatieve weerstandsverandering van $\frac{\Delta R}{R}\sim 10^{−6}$ worden waargenomen. Om deze gevoeligheid te kunnen benutten, moeten andere zaken die de weerstandsmeting kunnen beïnvloeden zo goed mogelijk buitenspel gezet worden.



## Voorbereiding

Lees de theorie en de methode sectie. Beantwoord de vragen.

<!-- #region -->
**Theorie:**  

Hier bespreken we eerst welke processen een rol spelen als we het hebben over de oploswarmte van zout in water. Theoretisch gezien is het oplossen is eigenlijk een twee-stapsproces. De eerste stap is het opbreken van de vaste stof, dat kost energie (energy of crystallization). In de tweede stap worden de deeltjes omgeven door watermoleculen. Deze stap levert energie op (energy of hydration). Bij het opbreken van stoffen zoals NaCl, zijn het de ionische verbindingen die worden verbroken. Deze energie is gelijk aan de energie die vrijkomt als de vaste stof uit de afzonderlijke ionen wordt gevormd. Vandaar de naam “energy of crystallization”.

```{figure} zout_ion.png
:label: fig_zout_ion
:width: 70%

Ion - dipool interactie tussen $\text{Na}^+$-ion en watermoleculen
```

Bij de omringing van een enkel ion door watermoleculen spelen ion-dipool krachten een rol. Ion-dipool krachten worden gevormd als de licht positieve waterstofatomen van het $\text{H}_2\text{O}$ molecuul worden aangetrokken door negatieve $\text{Cl}^-$ ionen en de licht negatieve zuurstofionen van het $\text{H}_2\text{O}$ molecuul worden aangetrokken door de positieve $\text{Na}^+$ ionen. Het vormen van deze omringing wordt “hydration” genoemd. De oploswarmte is het verschil van de kristallisatie-energie en de hydratatie-energie.

```{figure} zout_diagram.png
:label: fig_zout_diagram
:width: 70%

Energie (links) en entropie (rechts) schema voor de oplossen van NaCl in water.
```


In een energieschema zien de bijdragen er als volgt uit (@fig_zout_diagram). De energieverschillen zijn aangeduid met de grootheid $\Delta H$.
(s), (g) betekent ionen in resp. vaste stof en gas fase. (aq) staat voor ionen in water opgelost. We zien dat het oplossen van NaCl in water netto ongeveer $4 \mathrm{kJ/mol}$ kost. Het is het verschil van twee grote bijdragen. Het is een **endotherm proces**. Dit is zeker niet voor alle zouten het geval.

Het oplossen van zout in water verloopt spontaan als de vrije energie $G$ kleiner wordt. De verandering van de vrije energie (bijv. bij het oplossen) wordt gegeven door:

$$\Delta G = \Delta H - T\Delta S$$

$\Delta H$ is de verandering van de enthalpie (de totale warmte-inhoud van een systeem aangeeft bij een constante druk). In ons voorbeeld is dat energie die het kost om het zout op te lossen. $T$ is de absolute temperatuur. $\Delta S$ is de verandering van de entropie. We zien dat een vergroting van de wanorde, dus een toename van de entropie ($\Delta S > 0$) leidt tot een afname van $G$ ($\Delta G < 0$), de voorwaarde voor een spontaan proces. Hetzelfde geldt voor een afname van $H$. $\Delta H < 0$ betekent dat er warmte vrijkomt (exotherm). Uiteindelijk is het de balans die bepaalt of een proces spontaan verloopt of niet. Je kunt nu zelf inzien waarom het verhogen van de temperatuur vaak helpt.
<!-- #endregion -->

```{exercise}
:label: ex_zout_1

1. Laat zien dat de temperatuurverandering $\Delta T$ in Celsius van het water gegeven wordt door:
\begin{equation}
\Delta T = \Delta H^e \frac{m_{s}}{m_w C_p M_{s}}
\end{equation}

Waar $\Delta H^e$ de enthalpeverandering in kJ/mol is, $m_s$ de massa van het toegevoegde zout, $M_s$ de molaire massa van zout in kg/mol, $m_w$ de massa van water in kg, en $C_p$ de soortelijke warmte van water in J/kgK. Gebruik hiervoor de gegeven theorie.

2. Leid een uitdrukking af voor de onzekerheid $u(T)$ als we de onzekerheden in $C_p$ en $M_s$ buiten beschouwing laten omdat deze waarden nauwkeurig bekend zijn. De onzekerheden in $m_w$, $\Delta H^e$ en $m_s$ worden gegeven door $u(m_w)$, $u(\Delta H^e)$ en $u(m_s)$.

3. Bereken de temperatuursverandering van $100 \; \mathrm{mL}$ water als we daar $1 \; \mathrm{g}$ zout in oplossen.
```


```{solution} ex_zout_1

```

```python

```

**Meetinstrumenten en instellingen:**  

In de praktijk wordt veel gebruik gemaakt van weerstandsthermometers: thermometers waarvan de weerstandswaarde verandert als functie van de temperatuur. Hier maken we gebruik van een Pt100 sensor, een thermometer gemaakt van platinum. De voordelen van platina zijn ondermeer chemische stabiliteit, relatief gemakkelijk te maken, beschikbaarheid van de draad in zeer zuivere vorm en de zeer goede reproduceerbaarheid van de elektrische eigenschappen. De eigenschappen van een platina thermometer zijn vastgelegd in IEC and DIN specificaties. In @fig_zout_tempsensor is de weerstand van een Pt100 sensor geplot als functie van de temperatuur in $\degree \mathrm{C}$.

```{figure} zout_tempsensor.png
:label: fig_zout_tempsensor
:width: 70%

Weerstand van een Pt100 sensor als functie van de temperatuur in $\degree \mathrm{C}$.
```

In het gebied tussen 0 en $850 \degree \mathrm{C}$ wordt de weerstand van de thermometer gegeven door:

$$ R(T) = R_0(1 + AT + BT^2)$$

$R(T)$ is de weerstand bij temperatuur $T$ (in $\degree \mathrm{C}$). $R_0$ heet de nominale weerstand en is de weerstand bij $0 \degree \mathrm{C}$ en is in ons geval gelijk aan $100 \Omega$.
De coëfficiënten $A$ en $B$ worden gegeven door:\
$A = 3.9083 \times 10^{-3} \degree \mathrm{C^{-1}}$\
$B = -5.775 \times 10^{-7} \degree \mathrm{C^{-2}}$


```{exercise}
:label: ex_zout_2
1. Laat zien dat bij kamertemperatuur en een temperatuur verandering van $1 \; \degree\mathrm{C}$ verwaarlozing van de kwadratische term geen significante verandering in $R(T)$ veroorzaakt.

2. Laat zien dat ook voor een ijking van de thermometer op basis van ijswater en kokend water de verwaarlozing van de kwadratische term niet leidt tot een wezenlijke verandering in de waarde van $A$.

3. Schrijf een protocol (meetprocedure) voor het ijken van de thermometer en het vergelijken van de gevonden waarden.

4. Op basis van bovenstaande, schrijf een functie die de verandering in temperatuur uit rekent op basis van een begin en eindwaarde van de weerstand.
```


```{solution} ex_zout_2
```

```python
A = 3.9083e-3
B = -5.775e-7
R_nom = 100.0

def temperature(R):
    return (-A + np.sqrt(A**2 - 4*B*(1 - R/R_nom))) / (2*B)

def thermistor_temp(R_t_0, R_t_e):
    return temperature(R_t_e) - temperature(R_t_0)

print(temperature(100.0), temperature(138.505))
```

```{exercise}
:label: ex_zout_3
Schrijf op basis van de theorie en de `thermistor_temp` functie een functie die de oploswarmte van zout berekent.

```

```python
def oploswarmte(m_w, m_s, R_t_0, R_t_e):
    dT = thermistor_temp(R_t_0, R_t_e)
    return -(m_w/1000 * Cp * dT * Ms_eff) / (m_s/1000)
```

## Uitvoering

<!-- #region -->
**IJking**
Voer de bedachte ijking uit voor je thermometer.

**Methode**  
Voor het bepalen van een betrouwbare waarde van de oploswarmte van keukenzout doen we meerdere metingen steeds met een andere hoeveelheid water en zout. De hoeveelheid zout en water mag je zelf kiezen, maar bedenk dat de maximale concentratie $36 \; \mathrm{g}$ per $100 \; \mathrm{mL}$ (bij kamertemperatuur) is. Bedenk ook zelf of een grote hoeveelheid zout of een kleine hoeveelheid zout de voorkeur heeft. 

```{note} 
Zorg er wel voor dat de thermometer in een voldoende hoeveelheid water zit.
``` 

1. Meet de weerstandswaarde bij de gekozen hoeveelheid water. Let op! Na het inschenken van het water kan de weerstandswaarde fluctueren.
2. Voeg de gekozen waarde zout toe en meet de maximale weerstandsverandering.
3. Bepaal op basis van de weerstandsverandering de temperatuurverandering en daarmee de oploswarmte van zout met de bijbehorende onzekerheid.
4. Herhaal de procedure voor 5 verschillende hoeveelheden zout / water.


**Onzekerheidsbepaling**  
In het bepalen van een betrouwbare waarde van de oploswarmte van keukenzout maken we gebruik van een gewogen gemiddelde. Uit het [dictaat](https://contemporary-physicslab.github.io/FYPLC/notebook-6-measurement-uncertainty/#weighted-average) halen we dat het gewogen gemiddelde gegeven wordt door:
$$
    \overline{x} = \frac{w_1 x_1 + w_2 x_2 + w_3 x_3 + ... + w_n x_n}{w_1 + w_2 + w_3 + ... + w_n} = \frac{\sum_{i=1}^{N} w_i x_i}{\sum_{i=1}^{N} w_i}
$$

waarin $w_i$ een weegfactor is gegeven door $w_i=\frac{1}{u(x_i)^2}$. 

De onzekerheid in het gewogen gemiddelde wordt gegeven door:

$$ u(\overline{x}) = \sqrt{\frac{1}{\sum_i^N{w_i}}}$$




<!-- #endregion -->

**Onafhankelijke variable:** Massa van zout en water  
**Afhankelijke variable:**  Weerstand waarde
**Gecontroleerde variablen:**  Multimeter



**Setup (tekening of foto):**  



**Notities:**  





## Resultaten

**Metingen**

```python
m_s   = np.array([ 30.008,30.007 ,10.009 ,9.994 ,10.000]) # Massa van zout
m_w   = np.array([100 ,125.040,80.020 ,60.000 ,40.001 ]) # Massa van water
R_t_0 = np.array([108.522,108.170,108.228 ,108.250 ,108.555 ]) # Begin weerstand
R_t_e = np.array([116.028 ,114.265,111.407 ,112.525 ,114.670 ]) # Eind weerstand
```

```python
import numpy as np
import matplotlib.pyplot as plt

# constantes
A = 3.9083e-3
B = -5.775e-7
R_nom = 100.0
Cp = 4180          # J/(kg·K)

# calcium chloride 77% flakes (watervrij, maar 77% zuiver)
Ms_pure = 0.11098  # kg/mol (puur watervrij CaCl2)
purity = 0.77      # 77% zuiverheid
Ms_eff = Ms_pure / purity 

dT = thermistor_temp(R_t_0, R_t_e)
H = oploswarmte(m_w, m_s, R_t_0, R_t_e)

print("Temperatuurverandering (°C):", dT)
print("Oploswarmte (kJ/mol):", H/1000)

# onzekerheden
def u_R(R):
    return (0.00010*R + 0.00001*1000) / np.sqrt(3)

def u_temperature(R):
    T = temperature(R)
    return u_R(R) / abs(R_nom*(A + 2*B*T))

u_T_0 = u_temperature(R_t_0)
u_T_e = u_temperature(R_t_e)
u_dT = np.sqrt(u_T_0**2 + u_T_e**2)

# onzekerheid van de balans (B303-S) in kg
u_m = np.sqrt(0.001**2 + (0.001/np.sqrt(12))**2 + (0.002/np.sqrt(3))**2) / 1000

# onzekerheid in H 
u_H = np.abs(H) * np.sqrt(
    (u_dT/np.abs(dT))**2
    + (u_m/(m_s/1000))**2
    + (u_m/(m_w/1000))**2
)

# gewogen gemiddelde
w = 1/u_H**2
H_average = np.sum(w*H)/np.sum(w)
u_H_average = np.sqrt(1/np.sum(w))

print("Individuele onzekerheden (kJ/mol):", u_H/1000)
print("Gewogen gemiddelde (kJ/mol):", H_average/1000)
print("Onzekerheid (kJ/mol):", u_H_average/1000)
print(f"Eindresultaat: {H_average/1000:.2f} +/- {u_H_average/1000:.2f} kJ/mol")
```

**Notities:**  


**Verwerking**





```python jupyter={"outputs_hidden": true}

```

## Conclusie



## Andere notities, opmerkingen, verklaringen, ideeen etc
In dit experiment is de oploswarmte van calciumchloride (CaCl₂, 77% E509) bepaald op −36,0 ± 0,1 kJ/mol. Het negatieve teken toont aan dat het oplossen een exotherm proces is.

Deze waarde wijkt ca. 56% af van de literatuurwaarde van watervrij CaCl₂ (ca. −82,8 kJ/mol). Deze systematische fout wordt verklaard door:

Warmteverlies: De opstelling was niet geïsoleerd, waardoor warmte naar de omgeving weglekte en de gemeten temperatuurstijging te laag uitviel.

Onvolledige oplossing: Een deel van het zout was mogelijk nog niet opgelost op het moment dat de maximale temperatuur (via de weerstandssensor) werd bereikt.

Het experiment kan verbeterd worden door gebruik te maken van zuiverder zout (100%), een geïsoleerde calorimeter (dubbelwandige beker met deksel) en door minder zout te gebruiken waardoor het makkelijker en sneller oplost.



Klaar? Lever markdown file in.

```python

```

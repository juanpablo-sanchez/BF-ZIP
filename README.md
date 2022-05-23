# BF-ZIP
</head>

<body lang=ES style='tab-interval:35.4pt;word-wrap:break-word'>

<div class=WordSection1>

<p class=MsoNormal><span lang=ES-TRAD><o:p>&nbsp;</o:p></span></p>

<p class=MsoNormal><span lang=EN-US style='mso-ansi-language:EN-US'>PARAMETER
FILE FOR ZIP MODEL (ZIP_FB_h2_ratios_INV)<o:p></o:p></span></p>

<p class=MsoNormal style='margin-bottom:0cm'><span lang=EN-US style='mso-ansi-language:
EN-US'><o:p>&nbsp;</o:p></span></p>

<p class=MsoNormal style='margin-bottom:0cm'><span lang=EN-US style='mso-ansi-language:
EN-US'><o:p>&nbsp;</o:p></span></p>

<p class=MsoNormal style='margin-bottom:0cm'><span class=SpellE><span
lang=EN-US style='mso-ansi-language:EN-US'>d_rec</span></span><span lang=EN-US
style='mso-ansi-language:EN-US'><span style='mso-spacerun:yes'>              
</span># data file (trait, fixed factors, random factor, animal=<span
class=SpellE><span class=GramE>pedigree,Poisson</span>_average</span>)<o:p></o:p></span></p>

<p class=MsoNormal style='margin-bottom:0cm'><span class=SpellE><span
lang=EN-US style='mso-ansi-language:EN-US'>p_rec</span></span><span lang=EN-US
style='mso-ansi-language:EN-US'><span style='mso-spacerun:yes'>              
</span># pedigree file<o:p></o:p></span></p>

<p class=MsoNormal style='margin-bottom:0cm'><span lang=EN-US style='mso-ansi-language:
EN-US'>500<span style='mso-spacerun:yes'>                 </span># number of
rows in data file<o:p></o:p></span></p>

<p class=MsoNormal style='margin-bottom:0cm'><span lang=EN-US style='mso-ansi-language:
EN-US'>100<span style='mso-spacerun:yes'>                 </span># number of
rows in pedigree file <span class=SpellE>numero</span> de <span class=SpellE>animales</span>
<span class=SpellE>en</span> <span class=SpellE>fichero</span> pedigree<o:p></o:p></span></p>

<p class=MsoNormal style='margin-bottom:0cm'><span lang=EN-US style='mso-ansi-language:
EN-US'>5<span style='mso-spacerun:yes'>                   </span># number of
fixed and random (except the genetic <span class=SpellE>efect</span>) First
fixed them random<o:p></o:p></span></p>

<p class=MsoNormal style='margin-bottom:0cm'><span lang=EN-US style='mso-ansi-language:
EN-US'>5 6 2 196 189<span style='mso-spacerun:yes'>       </span># levels of
fixed and random factors<o:p></o:p></span></p>

<p class=MsoNormal style='margin-bottom:0cm'><span lang=EN-US style='mso-ansi-language:
EN-US'>0.0 0.0 0.0 0.1 0.1 # starting values for ratios of variance (when fixed
= 0)<o:p></o:p></span></p>

<p class=MsoNormal style='margin-bottom:0cm'><span lang=EN-US style='mso-ansi-language:
EN-US'>5.0<span style='mso-spacerun:yes'>                 </span># starting
value for Phenotypic variance (<span class=SpellE>log_ind_poisson_average</span>)<o:p></o:p></span></p>

<p class=MsoNormal style='margin-bottom:0cm'><span lang=EN-US style='mso-ansi-language:
EN-US'>0.10<span style='mso-spacerun:yes'>                </span># starting
value for heritability<o:p></o:p></span></p>

<p class=MsoNormal style='margin-bottom:0cm'><span lang=EN-US style='mso-ansi-language:
EN-US'>0.5<span style='mso-spacerun:yes'>                 </span># starting
value of the probability of Zero in the ZIP model <o:p></o:p></span></p>

<p class=MsoNormal style='margin-bottom:0cm'><span lang=EN-US style='mso-ansi-language:
EN-US'>2.0<span style='mso-spacerun:yes'>                 </span># Standard
deviation of the Metropolis Hasting steps for the update of <span class=SpellE>hte</span>
<span class=SpellE>log_lambda</span> parameters in the ZIP model <o:p></o:p></span></p>

<p class=MsoNormal style='margin-bottom:0cm'><span lang=EN-US style='mso-ansi-language:
EN-US'>10000<span style='mso-spacerun:yes'>               </span># number of
rounds<o:p></o:p></span></p>

<p class=MsoNormal style='margin-bottom:0cm'><span lang=EN-US style='mso-ansi-language:
EN-US'>1000<span style='mso-spacerun:yes'>                </span># burning<o:p></o:p></span></p>

<p class=MsoNormal style='margin-bottom:0cm'><span lang=EN-US style='mso-ansi-language:
EN-US'>1<span style='mso-spacerun:yes'>                   </span># number of
processors to use in the matrix multiplications and inversions<o:p></o:p></span></p>

<p class=MsoNormal style='margin-bottom:0cm'><span lang=EN-US style='mso-ansi-language:
EN-US'>0<span style='mso-spacerun:yes'>                   </span># set the <span
class=SpellE>heritabilityto</span> zero = 1, estimate the heritability = 0<o:p></o:p></span></p>

<p class=MsoNormal><span lang=EN-US style='mso-ansi-language:EN-US'><o:p>&nbsp;</o:p></span></p>

<p class=MsoNormal><span lang=EN-US style='mso-ansi-language:EN-US'><o:p>&nbsp;</o:p></span></p>

<p class=MsoNormal><span lang=EN-US style='mso-ansi-language:EN-US'>PARAMETER
FILE FOR NORMAL MODEL (FB_h2_ratios_INV)<o:p></o:p></span></p>

<p class=MsoNormal><span lang=EN-US style='mso-ansi-language:EN-US'><o:p>&nbsp;</o:p></span></p>

<p class=MsoNormal style='margin-bottom:0cm'><span lang=EN-US style='mso-ansi-language:
EN-US'><o:p>&nbsp;</o:p></span></p>

<p class=MsoNormal style='margin-bottom:0cm'><span class=SpellE><span
lang=EN-US style='mso-ansi-language:EN-US'>d_rec</span></span><span lang=EN-US
style='mso-ansi-language:EN-US'><span style='mso-spacerun:yes'>              
</span># data file (trait, fixed factors, random factor, animal=pedigree)<o:p></o:p></span></p>

<p class=MsoNormal style='margin-bottom:0cm'><span class=SpellE><span
lang=EN-US style='mso-ansi-language:EN-US'>p_rec</span></span><span lang=EN-US
style='mso-ansi-language:EN-US'><span style='mso-spacerun:yes'>              
</span># pedigree file<o:p></o:p></span></p>

<p class=MsoNormal style='margin-bottom:0cm'><span lang=EN-US style='mso-ansi-language:
EN-US'>500<span style='mso-spacerun:yes'>                 </span># number of
rows in data file<o:p></o:p></span></p>

<p class=MsoNormal style='margin-bottom:0cm'><span lang=EN-US style='mso-ansi-language:
EN-US'>100<span style='mso-spacerun:yes'>           </span><span
style='mso-spacerun:yes'>      </span># number of rows in pedigree file <span
class=SpellE>numero</span> de <span class=SpellE>animales</span> <span
class=SpellE>en</span> <span class=SpellE>fichero</span> pedigree<o:p></o:p></span></p>

<p class=MsoNormal style='margin-bottom:0cm'><span lang=EN-US style='mso-ansi-language:
EN-US'>5<span style='mso-spacerun:yes'>                   </span># number of
fixed and random (except the genetic <span class=SpellE>efect</span>) First
fixed them random<o:p></o:p></span></p>

<p class=MsoNormal style='margin-bottom:0cm'><span lang=EN-US style='mso-ansi-language:
EN-US'>5 6 2 196 189<span style='mso-spacerun:yes'>       </span># levels of
fixed and random factors<o:p></o:p></span></p>

<p class=MsoNormal style='margin-bottom:0cm'><span lang=EN-US style='mso-ansi-language:
EN-US'>0.0 0.0 0.0 0.1 0.1 # starting values for ratios of variance (when fixed
= 0)<o:p></o:p></span></p>

<p class=MsoNormal style='margin-bottom:0cm'><span lang=EN-US style='mso-ansi-language:
EN-US'>5.0<span style='mso-spacerun:yes'>                 </span># starting
value for Phenotypic variance <o:p></o:p></span></p>

<p class=MsoNormal style='margin-bottom:0cm'><span lang=EN-US style='mso-ansi-language:
EN-US'>0.10<span style='mso-spacerun:yes'>                </span># starting
value for heritability<o:p></o:p></span></p>

<p class=MsoNormal style='margin-bottom:0cm'><span lang=EN-US style='mso-ansi-language:
EN-US'>10000<span style='mso-spacerun:yes'>               </span># number of
rounds<o:p></o:p></span></p>

<p class=MsoNormal style='margin-bottom:0cm'><span lang=EN-US style='mso-ansi-language:
EN-US'>1000<span style='mso-spacerun:yes'>                </span># burning<o:p></o:p></span></p>

<p class=MsoNormal style='margin-bottom:0cm'><span lang=EN-US style='mso-ansi-language:
EN-US'>1<span style='mso-spacerun:yes'>                   </span># number of
processors to use in the matrix multiplications and inversions<o:p></o:p></span></p>

<p class=MsoNormal style='margin-bottom:0cm'><span lang=EN-US style='mso-ansi-language:
EN-US'>0<span style='mso-spacerun:yes'>                   </span># set the <span
class=SpellE>heritabilityto</span> zero = 1, estimate the heritability = 0<o:p></o:p></span></p>

<p class=MsoNormal><span lang=EN-US style='mso-ansi-language:EN-US'><o:p>&nbsp;</o:p></span></p>

<p class=MsoNormal><span lang=EN-US style='mso-ansi-language:EN-US'><o:p>&nbsp;</o:p></span></p>

<p class=MsoNormal><span lang=EN-US style='mso-ansi-language:EN-US'><o:p>&nbsp;</o:p></span></p>

<p class=MsoNormal><span lang=EN-US style='mso-ansi-language:EN-US'><o:p>&nbsp;</o:p></span></p>

</div>

</body>

</html>

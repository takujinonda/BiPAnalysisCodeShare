# To do

## DVL Steps

1. [x] Read in DVL tags
2. [x] Read in behaviour records (video)
3. [x] Calculate static/dynamic acceleration and pitch
4. [ ] Identify thresholds
   1. [ ] flight period estimation (prominence level?)
   2. [ ] pitch and ODBA variation
   3. [ ] repeat for $round(\frac{20}{24}n_{hours})$
5. [ ] Run through system with obtained thresholds
6. [ ] Test

## Algorithm steps

1. [x] Read in
2. [x] Speed/distance calculations
   1. [x] Remove erroneous GPS positions
   2. [x] Remove periods close to colony
   3. [x] For DVL only: remove minutes surround 'AT' labelled behaviours for spectrogram moving sum differences signal
3. [ ] Calculate flight periods
   1. [x] Generate a full acceleration feature calculator (saves memory?)
      1. [x] pass all acceleration columns
      2. [x] must understand which is longitudinal and which is dorsoventral
      3. [x] $acc_S$ & $acc_D$
      4. [x] Pitch and $p_{mn}$, moving mean over 10 seconds
      5. [x] ODBA and $ODBA_{mn}$, moving mean over 10 seconds
   2. [x] Extract frequency content intensity difference between $3 < x < 5$ and $x > 5$, henceforth $i_f$
   3. [ ] Ignore data where duration $< 120$ minutes
   4. [ ] Find indeces of 20 one minute periods where $i_f$ is largest minimum of 5 minutes apart
4. [ ] Extract typical flight thresholds
   1. [x] Find peaks / troughs and ensure signal starts with a peak (for difference calculation)
   2. [ ] Find overlap of flapping with estimated flight periods
   3. [ ] Calculate interpeak difference, $\Delta_z$
   4. [ ] Repeat interpeak difference calculation for pitch $\Delta_p$
   5. [ ] For all flight periods:
      1. [ ] Find $\frac{3}{2} g( \max \Delta_p)$ where $g(x)$ is the median
      2. [ ] $ g( \max{f(x_D)})$ where $f(x)$ is a moving variance over 2 seconds
      3. [ ] $ g( \max{ODBA_{mn}}) $
      4.  [ ] $ g( \min{p}) $
      5. [ ] $ g( \max{p_{mn}}) $
      6. [ ] $ g( \min{p_{mn}}) $
      7. [ ] $ g( h({p_{mn}})) $ where $h(x)$ is the variance of x
      8. [ ] $ g( \overline{p}) + 2 g( h(p))$ 
      9. [ ] $ g(\overline{p}) - 2g(h(p))$
      10. [ ] $ g(\min{p}) $
  1. [ ] Generate ethogram (same `length` as acceleration record)
     1. [ ] Set all as "Rest"
     2. [ ] Calculate moving mean of pitch and moving variance of $x_D$
     3. [ ] Identify $\Delta_z$ greater than threshold
     4. [ ] Group flaps within $0.5$s
     5. [ ] Group into flapping bouts (gaps $> 30$s)
     6. [ ] Set all flapping bouts as "Flight"
     7. [ ] Find large pitch changes $23.3$ seconds apart
     8. [ ] Where these align with "Flight" assignment and $ODBA_{mn} > ODBA_{fl}$, set as "Forage"
     9. [ ] Indicate dives
        1.  [ ] Search for preceeding dives
        2.  [ ] Check for large downward pitch change
        3.  [ ] Check if dive continues to reach above-normal pitch levels
   2. Set all foraging $<1$s or any $<2$s with two large downard pitch changes to `"Unknown"`


### Current System

* Read in and formatting
  * AxyTrek (`readAxyGPS`)
  * BiP-formatted (`readBiP`)
    * Format the BiP datetime (`dtFormat`)
* Distances and speed
  * Distances of multiple GPS points to single (`gps_distanceSingle`)
  * Distance and speed of GPS points and datetime (`gps_speed`)
  * Distance and speed of GPS points and datetime with speed threshold (`distSpeed`)
  * Remove acceleration and GPS data within threshold distance of capture site (`removeNearCS`)
* DSP
  * Generate and apply equiripple filter (`lowEquiFilt`)
  * Create a spectrogram with a 4 second hamming window with 85% overlap (`hammingSpect`)
  * Differences between summed spectrogram energies within defined frequencies (`rollingSpecSum`)
  * 
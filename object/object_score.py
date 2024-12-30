def compute_object_score(attributes):
    """
    To help compare multiple targets for astrophotography potential for a night or even period , we must calculate a
    score for each target and for each time segment we observe (see IntervalType.py for the options). If for example
    we are measuring the potential for a 5-hour astronomical night, and we fragment into half-hour intervals,
    the algorithm will compute 10 different scores for the target in that night and will add them up to create a
    master score for the night. When comparing to other targets' master score, we can easily decide which target is
    worth the effort and which is not.

    If a target is not visible or is below 25 degrees Altitude (too low for quality imagery), we break the
    computation with a 0 score.

    If the moon is visible but the conditions are not acceptable (moon phase too advanced or moon closer than 90
    degrees Azimuth to our target), then again we return a score of 0. Astr

    The closer the target is to Zenith at the observed time, the higher score we add up. If the target is within 10
    degreez Altitude of Zenith, we add 3 score points for exceptional imagery potential.

    Astrophotography done under an acceptable moon should also be rewarded, so if the moon is visible and the moon
    conditions are good, the target adds a higher score if its position is antipodal (other side of the sky dome) to
    the moon position. Narrowband signal acquisition is recommended during full moon to lower the effects of moon light.
    """

    # Mandatory filters.
    if not attributes["target_visible"]:
        return 0
    if not attributes["is_target_above_25_degree"]:
        return 0
    if attributes["is_moon_visible"] and not attributes["is_moon_phase_acceptable"]:
        return 0

    score = 1  # Default score if the target is visible and above 25 degrees.

    # Scoring based on target's position relative to zenith.
    if attributes["is_target_in_45_degrees_of_zenith"]:
        if attributes["is_target_in_10_degrees_of_zenith"]:
            score += 3
        elif attributes["is_target_in_30_degrees_of_zenith"]:
            score += 2
        else:
            score += 1

    # Scoring based on moon visibility and position.
    if attributes["is_moon_visible"]:
        if attributes["is_moon_phase_acceptable"] and not attributes["is_moon_in_180_degrees"]:
            score += 2
    else:
        score += 1

    return score

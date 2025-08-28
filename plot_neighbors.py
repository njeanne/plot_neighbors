#!/usr/bin/env python3

"""
Created on 26 Aug. 2025
"""

__author__ = "Nicolas JEANNE"
__copyright__ = "GNU General Public License"
__email__ = "jeanne.n@chu-toulouse.fr"
__version__ = "1.0.0"

import argparse
import logging
import os
import re
import statistics
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import yaml


def create_log(log_path, level, out_dir):
    """
    Create the log as a text file and as a stream.

    :param log_path: the path of the log.
    :type log_path: str
    :param level: the level og the log.
    :type level: str
    :param out_dir: the result directory path.
    :type out_dir: str
    """
    os.makedirs(out_dir, exist_ok=True)
    if not log_path:
        log_path = os.path.join(out_dir, f"{os.path.splitext(os.path.basename(__file__))[0]}.log")
    log_level_dict = {"DEBUG": logging.DEBUG,
                      "INFO": logging.INFO,
                      "WARNING": logging.WARNING,
                      "ERROR": logging.ERROR,
                      "CRITICAL": logging.CRITICAL}

    log_level = log_level_dict["INFO"]
    if level:
        log_level = log_level_dict[level]

    if os.path.exists(log_path):
        os.remove(log_path)

    logging.basicConfig(format="%(asctime)s %(levelname)s:\t%(message)s",
                        datefmt="%Y/%m/%d %H:%M:%S",
                        level=log_level,
                        handlers=[logging.FileHandler(log_path), logging.StreamHandler()])
    logging.info(f"version: {__version__}")
    logging.info(f"CMD: {' '.join(sys.argv)}")


def check_optional_args(args_domains, args_embedded_domains, args_residues_distance):
    """Check the combination of optional arguments

    :param args_domains: the domain's argument.
    :type args_domains: str
    :param args_embedded_domains: the embedded domain's argument.
    :type args_embedded_domains: bool
    :param args_residues_distance: the value of the residues distance argument.
    :type args_residues_distance: int
    """
    if args_embedded_domains and not args_domains:
        logging.warning("--embedded-domains will not be used as the --domains argument is missing.")
    if args_residues_distance != 4 and not args_domains:
        logging.warning(f"--residues-distance {args_residues_distance} will not be used as the --domains argument is "
                        f"missing.")


def get_domains(domains_path, use_embedded):
    """
    Load the domain file and fill in domains that are not covered.

    :param domains_path: the path to the domains file.
    :type domains_path: str
    :param use_embedded: if an embedded domain should be used as a domain or if only the embedding domain should be
    processed.
    :type use_embedded: bool
    :return: the filled-in domains data frame.
    :rtype: pd.Dataframe
    """
    logging.info(f"domains embedded in other domains will{' ' if use_embedded else ' not '}be used in the contacts "
                 f"by domain plot.")
    raw = pd.read_csv(domains_path, sep=",", header=0, names=["domain", "start", "stop", "color"], index_col=False)

    # check for embedded entries
    embedded_raw_idx = []
    previous_stop = 0
    for idx, row in raw.iterrows():
        if row["stop"] < previous_stop:
            embedded_raw_idx.append(idx)
        else:
            previous_stop = row["stop"]

    # update the domains
    data = {"domain": [], "start": [], "stop": [], "color": []}
    expected_start = 1
    previous = {"embedded": False, "domain": None, "stop": None, "color": None}
    for idx, row in raw.iterrows():
        if idx in embedded_raw_idx:  # case of embedded entry
            if use_embedded:
                # modify the previous end of the embedding domain
                data["stop"][-1] = row["start"] - 1
                # register the embedded domain
                data["domain"].append(row["domain"])
                data["start"].append(row["start"])
                data["stop"].append(row["stop"])
                data["color"].append(row["color"])
                # record the end of the domain where the embedded is
                previous["embedded"] = True
                expected_start = row["stop"] + 1
        else:
            if previous["embedded"]:
                # register the end of the domain where the previous was
                data["domain"].append(previous["domain"])
                data["start"].append(expected_start)
                data["stop"].append(previous["stop"])
                data["color"].append(previous["color"])
                expected_start = previous["stop"] + 1
                previous["embedded"] = False
            if row["start"] > expected_start:  # between domains
                # record the undefined domain
                if idx == 0:  # before first domain
                    data["domain"].append(f"before {row['domain']}")
                else:
                    data["domain"].append(f"between {previous['domain']} and {row['domain']}")
                data["start"].append(expected_start)
                data["stop"].append(row["start"] - 1)
                data["color"].append("#cecece")
            # record the current row domain
            data["domain"].append(row["domain"])
            data["start"].append(row["start"])
            data["stop"].append(row["stop"])
            data["color"].append(row["color"])
            previous["domain"] = row["domain"]
            previous["stop"] = row["stop"]
            previous["color"] = row["color"]
            expected_start = row["stop"] + 1

    # add a last row for the case of a contact after the last domain
    data["domain"].append(f"after {row['domain']}")
    data["start"].append(row["stop"] + 1)
    data["stop"].append(row["stop"] + 10000)
    data["color"].append(None)

    df = pd.DataFrame(data)
    return df


def extract_roi_id(domains, roi_coord):
    """
    Extract the Region Of Interest (ROI) identity from the domain file using the start and stop coordinates. If no
    match is performed, the ROI is returned as the start and stop coordinates.

    :param domains: the domain information.
    :type domains: pandas.dataframe
    :param roi_coord: the region of interest's start and stop coordinates.
    :type roi_coord: list
    :return: the ROI identity.
    :rtype: str
    """
    try:
        roi_id = domains.loc[(domains["start"] == roi_coord[0]) & (domains["stop"] == roi_coord[1])]["domain"].values[0]
    except IndexError:
        roi_id = f"{roi_coord[0]}-{roi_coord[1]}"
        logging.warning(f"no domains match with the coordinates {roi_id} in the domains CSV file provided, this "
                        f"coordinates will be used to named the Region Of Interest instead of a domain name.")
    return roi_id


def get_contacts_analysis_parameters(parameters_path):
    """
    Get the analysis parameters from the previous analysis from trajectories_contacts.py script.

    :param parameters_path: the path to the YAML parameters file.
    :type parameters_path: str
    :return: the parameters.
    :rtype: dict
    """
    with open(parameters_path, "r") as file_handler:
        parameters = yaml.safe_load(file_handler.read())
        logging.info("Parameters used for trajectory contacts search:")
        for p_key, p_value in parameters.items():
            if type(p_value) is dict:
                logging.info(f"\t{p_key}:")
                for p_key_2, p_value_2 in p_value.items():
                    logging.info(f"\t\t{p_key_2}:\t{p_value_2}")
            else:
                logging.info(f"\t{p_key}:\t{', '.join(p_value) if p_key == 'trajectory files processed' else p_value}")
    return parameters


def extract_roi(roi_to_extract):
    """
    Extract the region of interest (roi) start's and stop's coordinates.

    :param roi_to_extract: the coordinates ot the region of interest, as 100-200 i.e.,
    :type roi_to_extract: str
    :raises ArgumentTypeError: is not between 0.0 and 100.0
    :return: the region of interest (roi) start's and stop's coordinates.
    :rtype: list
    """
    pattern_roi_to_extract = re.compile("(\\d+)-(\\d+)")
    match_roi_to_extract = pattern_roi_to_extract.search(roi_to_extract)
    if match_roi_to_extract:
        roi_extracted = [int(match_roi_to_extract.group(1)), int(match_roi_to_extract.group(2))]
    else:
        raise argparse.ArgumentTypeError(f"'{roi_to_extract}' argument is malformed, it should be two integers "
                                         f"separated by an hyphen, i.e: '100-200'.")
    return roi_extracted


def outliers_neighbors(path_neighbors, proportion_thr, distance_thr):
    """
    Remove the neighbors contacts under the frames' proportions and the residues' distance thresholds.

    :param path_neighbors: the path to the neighbors file.
    :type path_neighbors: str
    :param proportion_thr: the neighbors contacts frames proportion's threshold.
    :type proportion_thr: float
    :param distance_thr: the residues' distance threshold.
    :type distance_thr: int
    :return: the dataframe of unique residues pairs contacts.
    :rtype: pd.Dataframe
    """
    df_init = pd.read_csv(path_neighbors, sep=",")
    # remove the neighborhood contacts under the frames' proportion threshold.
    df = df_init[df_init["proportion frames (%)"] >= proportion_thr]
    logging.debug(f"{len(df)}/{len(df_init)} neighborhood contacts present in {proportion_thr}% of the molecular "
                  f"dynamics frames.")
    # remove rows with too close distance between the residues
    idx_to_remove_for_residue_distance = []
    for idx, row in df.iterrows():
        if abs(row["residue 2 position"] - row["residue 1 position"]) < distance_thr:
            idx_to_remove_for_residue_distance.append(idx)
    df.drop(idx_to_remove_for_residue_distance, inplace=True, axis=0)
    # reset the index of the dataframe from 0
    df.reset_index(inplace=True, drop=True)
    logging.debug(f"{len(df)}/{len(df_init)} atoms contacts remaining with a minimal residues distance threshold of "
                  f"{distance_thr}.")

    return df


def update_domains(df, domains, out_dir, params, roi_id):
    """
    Get the domains for the acceptor and the donor in pairs.

    :param df: the dataframe of unique residues pairs contacts.
    :type df: pd.Dataframe
    :param domains: the domain's coordinates.
    :type domains: pd.Dataframe
    :param out_dir: the path output directory.
    :type out_dir: str
    :param params: the parameters used in the previous trajectory contacts analysis.
    :type params: dict
    :param roi_id: the region of interest name.
    :type roi_id: str
    :return: the pairs contacts dataframe updated with the regions.
    :rtype: pd.Dataframe
    """
    residue_1_regions = [None] * len(df)
    residue_2_regions = [None] * len(df)
    for idx, row_contacts in df.iterrows():
        for _, row_domains in domains.iterrows():
            if row_domains["start"] <= row_contacts["residue 1 position"] <= row_domains["stop"]:
                residue_1_regions[idx] = row_domains["domain"]
            if row_domains["start"] <= row_contacts["residue 2 position"] <= row_domains["stop"]:
                residue_2_regions[idx] = row_domains["domain"]
    df.insert(3, "residue 1 domain", pd.DataFrame(residue_1_regions))
    df.insert(6, "residue 2 domain", pd.DataFrame(residue_2_regions))
    out_path = os.path.join(out_dir, f"neighborhood_{params['sample'].replace(' ', '_')}_{roi_id}.csv")
    df.to_csv(out_path, index=False)
    logging.info(f"{roi_id} neighborhood's contacts saved: {out_path}")
    return df


def domains_involved(df, domains):
    """
    By domain, create the neighborhood contacts dataframes by atoms and by residues.

    :param df: the dataframe.
    :type df: pd.Dataframe
    :param domains: the domains.
    :type domains: pd.Dataframe
    :return: the domains dataframes of the contacts by atom and by residue.
    :rtype: pandas.Dataframe
    """
    # get the number of contacts by domain, then by residues' pair
    data = {}
    for _, row_domains in domains.iterrows():
        neighbors_in_domain = df[df["residue 2 domain"] == row_domains["domain"]]
        if len(neighbors_in_domain) > 0:
            data[row_domains["domain"]] = {}
            for _, row_neighbor in neighbors_in_domain.iterrows():
                residues_combination = f"{row_neighbor['residue 1 position']} - {row_neighbor['residue 2 position']}"
                if residues_combination not in data[row_domains["domain"]]:
                    data[row_domains["domain"]][residues_combination] = 0
                data[row_domains["domain"]][residues_combination] += 1
    # get the number of atoms in contact and the number of residues with atoms in contact
    data_by_atom = {}
    data_by_residue = {}
    for domain in data:
        data_by_residue[domain] = len(data[domain])
        nb_atoms = 0
        for residues_combination in data[domain]:
            nb_atoms += data[domain][residues_combination]
        data_by_atom[domain] = nb_atoms
    # save as dataframes
    source_by_atom = pd.DataFrame.from_dict({"domain": data_by_atom.keys(),
                                             "number of contacts": data_by_atom.values()})
    source_by_residue = pd.DataFrame.from_dict({"domain": data_by_residue.keys(),
                                                "number of contacts": data_by_residue.values()})

    return source_by_atom, source_by_residue


def plot_neighbors(source, out_dir, params, roi_id, fmt, res_dist, by_atom):
    """
    Create the neighborhood contacts plot by domains.

    :param source: the dataframe.
    :type source: pd.Dataframe
    :param out_dir: the path of the output directory.
    :type out_dir: str
    :param params: the parameters used in the previous trajectory contacts analysis.
    :type params: dict
    :param roi_id: the region of interest name.
    :type roi_id: str
    :param fmt: the format for the plot.
    :type fmt: str
    :param by_atom: the contacts are displayed by atoms.
    :type by_atom: bool
    :param res_dist: the maximal residues distance in the amino acids chain.
    :type res_dist: int
    """
    # set color and plot text values
    if by_atom:
        elt_type = "atom"
        plot_color="deeppink"
    else:
        elt_type = "residue"
        plot_color = "orangered"

    # set the seaborn plots style and size
    sns.set_style("darkgrid")
    sns.set_context("poster", rc={"grid.linewidth": 2})
    fig, ax = plt.subplots(figsize=(15, 15))
    sns.barplot(data=source, ax=ax, x="domain", y="number of contacts", color=plot_color)

    # modify the ticks labels for the X axis by adding new lines every 3 words
    modified_x_labels = [re.sub(r'(\w+ \w+ \w+)( )',
                                r'\1\n', x_label.get_text()) for x_label in ax.get_xticklabels()]
    # set the number of ticks for the X axis to avoid a matplotlib warning
    ax.set_xticks([num_tick for num_tick in range(len(modified_x_labels))])
    ax.set_xticklabels(modified_x_labels, rotation=45, horizontalalignment="right")

    ax.set_xlabel(None)

    ax.set_ylabel(f"{roi_id}: {elt_type} - {elt_type} neighborhood contacts", fontweight="bold")
    ax.text(x=0.5, y=1.1, s=f"{params['sample']}: {elt_type} - {elt_type} neighborhood contacts \nbetween {roi_id} and "
                            f"the protein domains",
            weight="bold", ha="center", va="bottom", transform=ax.transAxes)
    md_duration = f", MD: {params['parameters']['time']}" if "time" in params['parameters'] else ""
    ax.text(x=0.5, y=1.0,
            s=f"Maximal atoms distance: {params['parameters']['maximal atoms distance']} \u212B, minimal residues "
              f"distance: {res_dist}\n{params['parameters']['proportion contacts']}% of contacts in "
              f"{params['frames']} frames{md_duration}",
            alpha=0.75, ha="center", va="bottom", transform=ax.transAxes)
    path = os.path.join(out_dir,
                        f"neighborhood_{elt_type}_contacts_{params['sample'].replace(' ', '_')}_{roi_id}.{fmt}")
    fig.savefig(path, bbox_inches="tight")
    logging.info(f"Plot of {roi_id} neighborhood {elt_type} contacts by domain saved: {path}")


if __name__ == "__main__":
    descr = f"""
    {os.path.basename(__file__)} v. {__version__}

    Created by {__author__}.
    Contact: {__email__}
    {__copyright__}

    Distributed on an "AS IS" basis without warranties or conditions of any kind, either express or implied.

    From a CSV file of the neighborhood contacts during a molecular dynamics simulation and a YAML parameters file 
    produced by the script trajectories_neighbors.py (https://github.com/njeanne/trajectories_neighbors), a heatmap 
    representing the residues neighborhood contacts.
    
    A Region Of Interest (ROI) is defined with a range of amino acids selected in the protein, on the heatmap the 
    neighborhood contacts on the ordinate axis will be the ones belonging to this ROI.

    If a domains CSV file is used with the option "--domains", a plot and a CSV file of the neighborhood contacts by 
    domains will be produced. For this CSV, if some domains are embedded in other domains, they can be displayed in the 
    outputs with the option "--use-embedded".
    """
    parser = argparse.ArgumentParser(description=descr, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--out", required=True, type=str, help="the path to the output directory.")
    parser.add_argument("-d", "--domains", required=True, type=str, default="",
                        help="the path to the CSV domains file. A comma separated file, the first column is the "
                             "annotation name, the 2nd is the residue start coordinate, the 3rd is the residue end "
                             "coordinate, the last one is the color to apply in hexadecimal format. The coordinate are "
                             "1-indexed.")
    parser.add_argument("-p", "--parameters", required=True, type=str,
                        help="the path to the trajectory neighborhood contacts analysis parameters (the YAML file in "
                             "the results directory of the trajectory_neighbors.py script.")
    parser.add_argument("-r", "--roi", required=True, type=str,
                        help="the residue 1 Region Of Interest (ROI) amino acids coordinates, the format should be two "
                             "digits separated by an hyphen, i.e: '100-200'.")
    parser.add_argument("-f", "--format", required=False, default="svg",
                        choices=["eps", "jpg", "jpeg", "pdf", "pgf", "png", "ps", "raw", "svg", "svgz", "tif", "tiff"],
                        help="the output plots format: 'eps': 'Encapsulated Postscript', "
                             "'jpg': 'Joint Photographic Experts Group', 'jpeg': 'Joint Photographic Experts Group', "
                             "'pdf': 'Portable Document Format', 'pgf': 'PGF code for LaTeX', "
                             "'png': 'Portable Network Graphics', 'ps': 'Postscript', 'raw': 'Raw RGBA bitmap', "
                             "'rgba': 'Raw RGBA bitmap', 'svg': 'Scalable Vector Graphics', "
                             "'svgz': 'Scalable Vector Graphics', 'tif': 'Tagged Image File Format', "
                             "'tiff': 'Tagged Image File Format'. Default is 'svg'.")
    parser.add_argument("-x", "--residues-distance", required=False, type=int, default=4,
                        help="when 2 atoms of different residues are in contact, the minimal distance in number of "
                             "residues that should separate them for a long range interaction. Default is 4 residues, "
                             "the number of residues in an alpha helix.")
    parser.add_argument("-y", "--proportion", required=False, type=float, default=50.0,
                        help="the minimal percentage of frames where 2 atoms of different residues are neighbors, "
                             "default is 50%%.")
    parser.add_argument("-e", "--embedded-domains", required=False, action="store_true",
                        help="for the outliers plot of the neighborhood contacts between a specific domain and the "
                             "whole protein, use the domains embedded in another domain. In example, if the domain 2 "
                             "is in domain 1, the plot will represent the domain 1 as: domain-1 domain-2 domain-1. If "
                             "this option is not used only the domain 1 will be used in the plot.")
    parser.add_argument("-l", "--log", required=False, type=str,
                        help="the path for the log file. If this option is skipped, the log file is created in the "
                             "output directory.")
    parser.add_argument("--log-level", required=False, type=str,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="set the log level. If the option is skipped, log level is INFO.")
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument("input", type=str, help="the neighborhood contacts CSV file.")
    args = parser.parse_args()

    # create the logger
    create_log(args.log, args.log_level, args.out)
    check_optional_args(args.domains, args.embedded_domains, args.residues_distance)

    # get the Region Of Interest if specified
    roi_limits = extract_roi(args.roi)

    # load the contacts analysis parameters
    parameters_contacts_analysis = None
    try:
        parameters_contacts_analysis = get_contacts_analysis_parameters(args.parameters)
    except ImportError as exc:
        logging.error(exc, exc_info=True)
        sys.exit(1)

    # load and format the domains' file
    domains_data = get_domains(args.domains, args.embedded_domains)
    # match the Region Of Interest coordinates with a domain
    region_of_interest = extract_roi_id(domains_data, roi_limits)

    outliers = outliers_neighbors(args.input, args.proportion, args.residues_distance)
    logging.info(f"{len(outliers)} unique residues pairs contacts (<= "
                 f"{parameters_contacts_analysis['parameters']['maximal atoms distance']} \u212B) with a distance of "
                 f"at least {args.residues_distance} residues"
                 f"{' in the region of interest '+args.roi if args.roi else ''} (residues pair may have multiple atoms "
                 f"contacts).")

    # get the neighborhood contacts
    neighborhood_contacts = update_domains(outliers, domains_data, args.out, parameters_contacts_analysis,
                                           region_of_interest)
    # get the neighborhood contacts by atoms and by residues
    by_atom, by_residue = domains_involved(neighborhood_contacts, domains_data)

    # plot neighborhood contacts by atom
    plot_neighbors(by_atom, args.out, parameters_contacts_analysis, region_of_interest, args.format,
                   args.residues_distance, by_atom=True)

    # plot neighborhood contacts by residue
    plot_neighbors(by_residue, args.out, parameters_contacts_analysis, region_of_interest, args.format,
                   args.residues_distance, by_atom=False)

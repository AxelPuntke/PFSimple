#include "PFSimpleTask.hpp"

#include "ConverterIn.hpp"
#include "ConverterOut.hpp"

#include "AnalysisTree/PlainTreeFiller.hpp"
#include "AnalysisTree/TaskManager.hpp"

#include <iostream>
#include <getopt.h>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <nlohmann/json.hpp>
#include <regex>

using namespace AnalysisTree;
using json = nlohmann::json;

#define LIST_CONTAINS(list, value) std::find(list.begin(), list.end(), value) != list.end()
void set_json_value(json& j, const std::string& path, const std::string& value);

nlohmann::json load_config(const std::string& filepath) {
    std::ifstream file_stream(filepath);
    if (!file_stream.is_open()) {
        throw std::runtime_error("Failed to open file: " + filepath);
    }

    try {
        return nlohmann::json::parse(file_stream, nullptr, true, true); // last 'true' enables comment support
    } catch (const nlohmann::json::parse_error& e) {
        throw std::runtime_error("JSON parse error in file " + filepath + ": " + e.what());
    }
}

int main(int argc, char** argv) {
    std::string input_file;
    std::string config_file;
    std::string output_file = "PFSimpleOutput.root";
    std::string plain_output_file = "PFSimplePlainTree.root";
    std::vector<std::pair<std::string, std::string>> overrides;

    //Parse command line args
    static struct option long_options[] = {
        {"output",       required_argument, nullptr, 'o'},
        {"plain-output", required_argument, nullptr, 'p'},
        {"set",          required_argument, nullptr, 's'},
        {0, 0, 0, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "o:p:s:", long_options, nullptr)) != -1) {
        switch (opt) {
            case 'o': output_file = optarg; break;
            case 'p': plain_output_file = optarg; break;
            case 's': {
                std::string arg = optarg;
                size_t eqPos = arg.find('=');
                if (eqPos != std::string::npos) {
                    std::string key = arg.substr(0, eqPos);
                    std::string value = arg.substr(eqPos + 1);
                    overrides.emplace_back(key, value);  // Save for later
                } else {
                    std::cerr << "Invalid --set format: use key=value\n";
                    return 1;
                }
                break;
            }
            default:
                return 1;
        }
    }

    if (optind + 2 > argc) {
        std::cerr << "Usage: " << argv[0] << " filelist.txt config.json [--output file] [--plain-output file] [--set key=value]\n";
        return 1;
    }

    input_file = argv[optind];
    config_file = argv[optind + 1];

    // Load config JSON file
    json config = load_config(config_file);

    // Apply deferred overrides
    for (const auto& [key, value] : overrides) {
        set_json_value(config, key, value);
    }

    // For demonstration
    //std::cout << "Modified config:\n" << config.dump(2) << std::endl;

    // Setup PFSimple from JSON config
    auto* man = TaskManager::GetInstance();

    auto* in_converter = new ConverterIn();
    in_converter->SetRecEventHeaderName(config["io"].value("receventheader_branchname", "RecEventHeader"));
    in_converter->SetRecTracksName(config["io"]["rectracks_branchname"]);
    in_converter->SetSimTracksName("SimParticles");
    in_converter->SetTrackCuts(new Cuts("Cut to reproduce KFPF", {EqualsCut((config["io"]["rectracks_branchname"].get<std::string>() + ".pass_cuts").c_str(), 1)}));
    in_converter->SetMcInfoAvailable(config["io"].value("mc_info_available", true));
    // Optional MC-match weight threshold for pid_mode 1 (0 = off, standard).
    // Set to 0.7 in the AFKF comparison config to replicate AFKF's track selection.
    in_converter->SetMcMatchWeightMin(config["io"].value("mc_match_weight_min", 0.0));

    std::vector<Decay> decays;

    for (const auto& decay_cfg : config["decays"]) {
        printf("Configuring decay %s\n", decay_cfg["mother"]["name"].get<std::string>().c_str());

        // Load daughters & cuts
        auto daughters_cfg = decay_cfg["daughters"];
        std::vector<Daughter> daughters;
        for (const auto& daughter_cfg : daughters_cfg) {
            std::vector<Pdg_t> pdg_codes = daughter_cfg["pdg_code"].get<std::vector<Pdg_t>>();

            if (pdg_codes.size() > 1)
                daughters.emplace_back(Daughter(pdg_codes[0], pdg_codes));
            else
                daughters.emplace_back(Daughter(pdg_codes[0]));

            auto& d = daughters.back();
            d.CancelCuts();
            if (daughter_cfg.contains("cuts")) {
              auto cuts = daughter_cfg["cuts"];
              if (cuts.contains("chi2prim")) d.SetCutChi2Prim(cuts["chi2prim"]);
              if (cuts.contains("cos"))      d.SetCutCos(cuts["cos"]);
            }
        }

        // Load mother
        auto mother_cfg = decay_cfg["mother"];
        Mother mother(mother_cfg["pdg_code"]);
        mother.CancelCuts();
        if (mother_cfg.contains("mass"))       mother.SetMassPdg(mother_cfg["mass"]);
        if (mother_cfg.contains("mass_sigma")) mother.SetMassPdgSigma(mother_cfg["mass_sigma"]);

        if (mother_cfg.contains("cuts")) {
            auto cuts = mother_cfg["cuts"];
            if (cuts.contains("dist"))       mother.SetCutDistance(cuts["dist"]);
            if (cuts.contains("distSV"))     mother.SetCutDistanceToSV(cuts["distSV"]);
            if (cuts.contains("chi2geo"))    mother.SetCutChi2Geo(cuts["chi2geo"]);
            if (cuts.contains("cosopen"))    mother.SetCutCosOpen(cuts["cosopen"]);
            if (cuts.contains("chi2topo"))   mother.SetCutChi2Topo(cuts["chi2topo"]);
            if (cuts.contains("costopo"))    mother.SetCutCosTopo(cuts["costopo"]);
            if (cuts.contains("LdL"))        mother.SetCutLdL(cuts["LdL"]);
            if (cuts.contains("L"))          mother.SetCutDecayLength(cuts["L"]);
            if (cuts.contains("distPVline")) mother.SetCutDistancePVLine(cuts["distPVline"]);
            if (cuts.contains("invmass"))    mother.SetCutInvMass(cuts["invmass"]);
        }

        // Secondary mother cuts (SM)
        if (decay_cfg.contains("secondary_mother_cuts")) {
            auto combo_cuts = decay_cfg["secondary_mother_cuts"];

            if (combo_cuts.contains("chi2geo"))  mother.SetCutChi2GeoSM(combo_cuts["chi2geo"].get<std::vector<Float_t>>());
            if (combo_cuts.contains("cosopen"))  mother.SetCutCosOpenSM(combo_cuts["cosopen"].get<std::vector<Float_t>>());
            if (combo_cuts.contains("chi2topo")) mother.SetCutChi2TopoSM(combo_cuts["chi2topo"].get<std::vector<Float_t>>());
            if (combo_cuts.contains("costopo"))  mother.SetCutCosTopoSM(combo_cuts["costopo"].get<std::vector<Float_t>>());
        }

        // Construct decay and set save options
        Decay decay(mother_cfg["name"], mother, {daughters});

        if (mother_cfg.contains("save_options")) {
            const auto& options = mother_cfg["save_options"];
            if (LIST_CONTAINS(options, "mass_constraint"))    decay.SetIsApplyMassConstraint();
            if (LIST_CONTAINS(options, "transport_to_pv"))    decay.SetIsTransportToPV();
            if (LIST_CONTAINS(options, "do_not_save_mother")) decay.SetIsDoNotWriteMother();
        }

        decays.push_back(decay);
    }

    // PID and IO settings
    auto pid_mode = config["pid_mode"];
    auto pid_purity_cfg = config["pid_purity"];

    in_converter->SetPidMode(pid_mode);
    if (pid_purity_cfg.contains("default"))    in_converter->SetPidPurity(pid_purity_cfg["default"]);
    if (pid_purity_cfg.contains("protons"))    in_converter->SetPidPurityProton(pid_purity_cfg["protons"]);
    if (pid_purity_cfg.contains("pions"))      in_converter->SetPidPurityPion(pid_purity_cfg["pions"]);
    if (pid_purity_cfg.contains("kaons"))      in_converter->SetPidPurityKaon(pid_purity_cfg["kaons"]);
    if (pid_purity_cfg.contains("deuterons"))  in_converter->SetPidPurityDeuteron(pid_purity_cfg["deuterons"]);
    if (pid_purity_cfg.contains("background")) in_converter->SetPidPurityBG(pid_purity_cfg["background"]);

    auto* pf_task = new PFSimpleTask();
    pf_task->SetInTask(in_converter);
    pf_task->SetDecays({decays});

    auto* out_converter = new ConverterOut();
    man->SetOutputName(output_file, "pTree");
    out_converter->SetSimEventHeaderName("SimEventHeader");
    out_converter->SetRecTracksName(config["io"]["rectracks_branchname"]);
    out_converter->SetSimTracksName("SimParticles");
    out_converter->SetPFSimpleTask(pf_task);
    out_converter->SetDecays(decays);
    if(LIST_CONTAINS(config["io"]["save_options"], "write_detailed_bg")) out_converter->SetIsWriteDetailedBG(true);
    out_converter->SetMcInfoAvailable(config["io"].value("mc_info_available", true));

    // Event-mixing background estimation (optional, top-level config block).
    if (config.contains("event_mixing")) {
        const auto& em_cfg = config["event_mixing"];
        const bool em_enabled = em_cfg.value("enabled", true);
        const int em_mix_pdg = em_cfg.value("mix_pdg", -211);
        const int em_buffer_size = em_cfg.value("mix_buffer_size", 5);
        pf_task->SetIsEventMixing(em_enabled);
        pf_task->SetMixPdg(em_mix_pdg);
        pf_task->SetMixBufferSize(em_buffer_size);
        out_converter->SetIsWriteEventMixing(em_enabled);
        if (em_cfg.contains("branchname"))
            out_converter->SetMixedCandidatesBranchName(em_cfg["branchname"].get<std::string>());

        if (em_enabled) {
            // Print a one-time summary of how event mixing applies to each configured
            // mother. Mixing is track-level and single-species: only daughters whose
            // PDG *exactly* matches `mix_pdg` (sign-sensitive, mirroring AFKF's
            // exact-PDG mix pools) are replaced by tracks from other events. The
            // current event is always excluded from its own buffer (read-then-fill),
            // so its real daughters never combine with each other.
            //
            //   - all daughters match  -> all mixed species  (every daughter is swapped)
            //   - some daughters match  -> partial mixing  (only the matching species swapped)
            //   - no daughter matches   -> not mixed  (event-mixing settings do not apply)
            //
            // NOTE: "all mixed species" does NOT guarantee every daughter comes from a
            // *different* event. The buffer holds up to `mix_buffer_size` previous
            // events, and their mixed-species tracks are combined into one flat pool,
            // so a combinatorially paired set of same-species daughters can occasionally
            // both come from the *same* buffered event (identical to AFKF). This leaves
            // a small residual of cross-event real-decay leakage, but no current-event
            // signal.
            fprintf(stderr, "\n=== Event mixing summary (mix_pdg=%d, buffer=%d, 9 centrality bins) ===\n",
                   em_mix_pdg, em_buffer_size);
            for (const auto& decay : decays) {
                const int n_daughters = decay.GetNDaughters();
                int n_match = 0;
                for (const auto& daughter : decay.GetDaughters()) {
                    const auto& pids = daughter.GetPids();
                    if (std::find(pids.begin(), pids.end(), em_mix_pdg) != pids.end())
                        ++n_match;
                }
                const char* status;
                if (n_match == 0) {
                    status = "NOT mixed";
                    fprintf(stderr, "WARNING: event-mixing settings do not apply to mother '%s' "
                            "(no daughter PDG matches mix_pdg=%d; no mixed background will be produced for it).\n",
                            decay.GetName().c_str(), em_mix_pdg);
                } else if (n_match < n_daughters) {
                    status = "PARTIAL";
                    fprintf(stderr, "WARNING: mother '%s' mixed background is only a partially mixed sample "
                            "(%d/%d daughters match mix_pdg=%d). A fully mixed implementation is not yet ready.\n",
                            decay.GetName().c_str(), n_match, n_daughters, em_mix_pdg);
                } else {
                    status = "ALL species mixed";
                }
                fprintf(stderr, "  %-20s : %s (%d/%d daughters match mix_pdg=%d)\n",
                       decay.GetName().c_str(), status, n_match, n_daughters, em_mix_pdg);
            }
            fprintf(stderr, "  note: mixed daughters are drawn from the previous %d events in the same\n"
                   "  centrality bin (flat pool); same-species daughters can occasionally both\n"
                   "  come from one buffered event (as in AFKF) - a small cross-event leak, not\n"
                   "  a current-event signal. The current event is excluded from its own buffer.\n",
                   em_buffer_size);
            fprintf(stderr, "===\n\n");
        }
    }

    std::vector<AnalysisTree::SimpleCut> vec_output_cuts = {};
    if (config.contains("output_cuts"))
    {
        for (const auto& cut_cfg : config["output_cuts"]) {
          if (cut_cfg.contains("equal"))
            vec_output_cuts.push_back(AnalysisTree::EqualsCut("Candidates." + cut_cfg["var"].get<std::string>(), cut_cfg["equal"]));
          else if (cut_cfg.contains("from"))
            vec_output_cuts.push_back(AnalysisTree::RangeCut("Candidates." + cut_cfg["var"].get<std::string>(), cut_cfg["from"], cut_cfg["to"]));
        }
    }

    if(LIST_CONTAINS(config["io"]["save_options"], "signal_only"))
        vec_output_cuts.push_back(AnalysisTree::SimpleCut({"Candidates.generation"}, []( std::vector<double>& var ) { return var.at(0) != 0; }));

    if (vec_output_cuts.size() > 0)
    {
        AnalysisTree::Cuts* output_cuts = new AnalysisTree::Cuts("output_cuts", vec_output_cuts);
        out_converter->SetOutputCuts(output_cuts);
    }
    
    man->AddTask(in_converter);
    man->AddTask(pf_task);
    man->AddTask(out_converter);
    man->SetVerbosityPeriod(100);
    
    man->Init({input_file}, {config["io"]["input_treename"]});
    man->Run(config["io"]["n_events"]);
    man->Finish();
    man->ClearTasks();

    if(LIST_CONTAINS(config["io"]["save_options"], "make_plain_tree")) {
        std::ofstream filelist;
        filelist.open("filelist.txt");
        filelist << output_file;
        filelist << "\n";
        filelist.close();

        auto* tree_task = new PlainTreeFiller();
        tree_task->SetOutputName(plain_output_file, "plain_tree");
        std::string branchname_rec = "Candidates";
        tree_task->SetInputBranchNames({branchname_rec});
        tree_task->AddBranch(branchname_rec);
        //tree_task->SetIsPrependLeavesWithBranchName(false); //Uncomment when most recent AnalysisTree is installed in cbmroot

        man->AddTask(tree_task);

        man->Init({"filelist.txt"}, {"pTree"});
        man->Run(config["io"]["n_events"]);
        man->Finish();
    }

    return 0;
}


/* 
Function to overide a value inside the json config, e.g. decays[0].mother.cuts.LdL=4.0
Command line argument usage example for lists: --set 'io.save_options=["signal_only"]'
*/
void set_json_value(json& config, const std::string& path, const std::string& value) {
    std::cout << "[DEBUG] Setting " << path << " to " << value << std::endl;

    std::regex token_regex(R"(([^.\[\]]+)(?:\[(\d+)\])?)");
    auto tokens_begin = std::sregex_iterator(path.begin(), path.end(), token_regex);
    auto tokens_end = std::sregex_iterator();

    json* current = &config;

    for (auto it = tokens_begin; it != tokens_end; ++it) {
        std::string key = (*it)[1];
        std::string index_str = (*it)[2];
        bool is_last = (std::next(it) == tokens_end);

        if (!index_str.empty()) {
            // Handle arrays
            int index = std::stoi(index_str);
            if (!(*current)[key].is_array()) {
                (*current)[key] = json::array();
            }

            while ((*current)[key].size() <= index) {
                (*current)[key].push_back(nullptr);
            }

            if (is_last) {
                json& target = (*current)[key][index];
                try {
                    if (value == "true" || value == "false") {
                        target = (value == "true");
                    } else if (value == "null") {
                        target = json();
                    } else if (!value.empty() && value.front() == '[' && value.back() == ']') {
                        target = json::parse(value);  // Parse list
                    } else if (value.find('.') != std::string::npos) {
                        target = std::stod(value);
                    } else {
                        target = std::stoi(value);
                    }
                } catch (...) {
                    target = value;
                }
                return;
            } else {
                current = &((*current)[key][index]);
            }
        } else {
            if (is_last) {
                json& target = (*current)[key];
                try {
                    if (value == "true" || value == "false") {
                        target = (value == "true");
                    } else if (value == "null") {
                        target = json();
                    } else if (!value.empty() && value.front() == '[' && value.back() == ']') {
                        target = json::parse(value);  // Parse list
                    } else if (value.find('.') != std::string::npos) {
                        target = std::stod(value);
                    } else {
                        target = std::stoi(value);
                    }
                } catch (...) {
                    target = value;
                }
                return;
            } else {
                current = &((*current)[key]);
            }
        }
    }
}
(() => {
  const model = "Qwen2.5-1.5B-Instruct-q4f16_1-MLC";
  const webllmUrl = "https://cdn.jsdelivr.net/npm/@mlc-ai/web-llm@0.2.89/+esm";
  const button = document.getElementById("explain_trial_qwen");
  const status = document.getElementById("trial_qwen_status");
  const answer = document.getElementById("trial_qwen_answer");
  const progressPanel = document.getElementById("trial_qwen_progress_panel");
  const progressBar = document.getElementById("trial_qwen_progress");
  const progressLabel = document.getElementById("trial_qwen_progress_label");
  let engine;
  let busy = false;
  const regimen = document.getElementById("combo");
  const evidence = JSON.parse(document.getElementById("trial_qwen_facts").dataset.facts);
  let selected = regimen.value;

  function clearOldExplanation() {
    if (selected === regimen.value) return;
    selected = regimen.value;
    answer.textContent = "";
    progressPanel.hidden = true;
    status.textContent = busy ? "Finishing the previous request." : "";
  }

  function showRegistryFacts(facts) {
    const counts = document.createElement("div");
    counts.textContent =
      `${facts.regimen} appears in ${facts.arms} matched arms across ${facts.studies} studies. ` +
      `Radiation was listed in ${facts.radiation_arms} arms, additional drugs or biological agents in ` +
      `${facts.additional_agent_arms}, and other interventions in ${facts.other_intervention_arms}. Counts may overlap.`;
    answer.append(counts);

    const examples = document.createElement("ul");
    for (const example of facts.examples || []) {
      if (!/^NCT[0-9]{8}$/.test(example.NCTId)) continue;
      const item = document.createElement("li");
      const link = document.createElement("a");
      link.href = `https://clinicaltrials.gov/study/${example.NCTId}`;
      link.textContent = example.NCTId;
      item.append(link, ` — ${example.Arm}: ${example.AdditionalEntries}`);
      examples.append(item);
    }
    if (examples.childElementCount) answer.append(examples);
  }

  // Selectize emits a jQuery change on the original select.
  regimen.addEventListener("change", clearOldExplanation);
  window.jQuery?.(regimen).on("change", clearOldExplanation);
  button.disabled = false;
  status.textContent = "";

  button.addEventListener("click", async () => {
    if (busy) return;
    clearOldExplanation();
    const current = regimen.value;
    const facts = evidence[current];
    answer.textContent = "";
    progressPanel.hidden = true;
    if (!facts?.arms) {
      status.textContent = "No matched registry arms are available for this regimen.";
      return;
    }
    showRegistryFacts(facts);
    if (!navigator.gpu) {
      status.textContent = "Qwen requires a WebGPU-compatible browser. Use Summarise registry evidence for the local summary.";
      return;
    }

    busy = true;
    button.disabled = true;
    try {
      if (!engine) {
        progressPanel.hidden = false;
        progressBar.removeAttribute("value");
        progressLabel.textContent = "Starting…";
        status.textContent = "Preparing the browser runtime.";
        const { CreateMLCEngine } = await import(webllmUrl);
        if (current === regimen.value) status.textContent = "Downloading or loading model files.";
        engine = await CreateMLCEngine(model, {
          initProgressCallback: ({ progress }) => {
            if (current !== regimen.value || !Number.isFinite(progress)) return;
            const percent = Math.max(0, Math.min(100, Math.floor(progress * 100)));
            progressBar.value = percent;
            progressLabel.textContent = `${percent}%`;
            if (percent === 100) status.textContent = "Initialising the model.";
          },
        });
      }
      if (current !== regimen.value) return;
      if (!progressPanel.hidden) {
        progressBar.value = 100;
        progressLabel.textContent = "100% · Model ready";
      }

      status.textContent = "Writing the explanation…";
      const response = await engine.chat.completions.create({
        messages: [
          { role: "system", content: [
            "Write two short sentences interpreting what matched cervical cancer trial registry arms show.",
            "The verified counts and links are already displayed. Do not write counts, trial IDs,",
            "drug names, or facts about individual studies.",
            "Explain that additional interventions can complicate interpretation of the regimen alone",
            "when present. Registry matches show what was studied and cannot establish efficacy.",
            "Write plain text without treatment recommendations.",
          ].join(" ") },
          { role: "user", content: JSON.stringify({
            hasRadiation: facts.radiation_arms > 0,
            hasAdditionalAgents: facts.additional_agent_arms > 0,
            hasOtherInterventions: facts.other_intervention_arms > 0,
          }) },
        ],
        temperature: 0,
        max_tokens: 280,
      });
      if (current !== regimen.value) return;

      const text = response.choices[0]?.message?.content?.trim() || "";
      if (!text || /\bNCT\s*\d+|[0-9]/i.test(text)) {
        status.textContent = "The verified registry summary is shown below. Qwen did not return a usable interpretation.";
        return;
      }
      const draft = document.createElement("div");
      draft.textContent = text;
      answer.prepend(draft);
      status.textContent = "Qwen interpretation and verified registry references";
    } catch {
      if (current === regimen.value) {
        progressPanel.hidden = true;
        status.textContent = "Qwen could not run on this device. Use Summarise registry evidence for the local summary.";
      }
    } finally {
      busy = false;
      button.disabled = false;
      if (current !== regimen.value) status.textContent = "";
    }
  });
})();

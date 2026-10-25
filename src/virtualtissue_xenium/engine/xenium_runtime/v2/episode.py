"""Validate immutable recordings and produce a portable, inference-free player."""
import base64
import gzip
import hashlib
from pathlib import Path

from ..io import canonical,fingerprint


def engine_fingerprint():
    """Bind a paid journal to the exact implementation, not only its parameters."""
    source=Path(__file__).resolve().parents[1]
    return fingerprint({p.relative_to(source).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(source.rglob("*.py"))})


def validate_episode(episode,atlas=None):
    if not isinstance(episode,dict) or episode.get("format") != "luad-episode-v2":
        raise ValueError("Expected an autonomous v2 recording")
    body = {k:v for k,v in episode.items() if k != "fingerprint"}
    if fingerprint(body) != episode.get("fingerprint"):
        raise ValueError("Recording fingerprint mismatch")
    if fingerprint(episode["final_state"]) != episode["final_state_fingerprint"]:
        raise ValueError("Final state fingerprint mismatch")
    if atlas is not None and atlas["fingerprint"] != episode["atlas_fingerprint"]:
        raise ValueError("Recording and measured atlas differ")
    frames = episode.get("frames",[])
    if not frames or frames[-1]["time_min"] != episode["summary"]["time_min"]:
        raise ValueError("Final frame missing")
    final = episode["final_state"]
    present = [i for i,v in enumerate(final["present"]) if v]
    last = frames[-1]
    if last["ids"] != [final["ids"][i] for i in present]:
        raise ValueError("Final frame identities mismatch")
    for j,i in enumerate(present):
        if [last["x"][j],last["y"][j]] != final["xy"][i]:
            raise ValueError("Final frame positions mismatch")
        if last["states"][j]["viability"] != final["state"][i]["viability"]:
            raise ValueError("Final frame viability mismatch")
    return dict(verified=True,frames=len(frames),final_state_fingerprint=episode["final_state_fingerprint"])


def portable_player(episode,atlas,web_dir,destination):
    validate_episode(episode,atlas)
    web_dir = Path(web_dir)
    raw = gzip.compress(canonical(dict(episode=episode,atlas=atlas)),mtime=0)
    data = base64.b64encode(raw).decode("ascii")
    digest = hashlib.sha256(raw).hexdigest()
    html = (web_dir/"autonomous.html").read_text(encoding="utf8")
    html = html.replace('<link rel="stylesheet" href="/autonomous.css">',"<style>"+(web_dir/"autonomous.css").read_text(encoding="utf8")+"</style>")
    script = (web_dir/"autonomous.js").read_text(encoding="utf8")
    html = html.replace('<script src="/autonomous.js" defer></script>',f'<script id="portable-data" type="application/octet-stream" data-sha256="{digest}">{data}</script>')
    html = html.replace("</body>","<script>"+script+"</script></body>")
    Path(destination).write_text(html,encoding="utf8")

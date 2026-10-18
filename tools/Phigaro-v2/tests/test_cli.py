import pytest

from phigaro import bounded


PREDICT_ARGS = [
    "predict",
    "--input",
    "input.fa",
    "--output",
    "output",
    "--models-hmm",
    "models.hmm",
    "--models-sha256",
    "0" * 64,
]


def test_predict_help_keeps_timeout_concise(capsys):
    with pytest.raises(SystemExit) as exc:
        bounded.main(["predict", "--help"])
    assert exc.value.code == 0
    help_text = capsys.readouterr().out
    assert len(help_text) < 2000
    assert "--timeout-seconds SECONDS" in help_text
    assert "1-14400; default: 3600" in help_text


@pytest.mark.parametrize(
    "timeout,expected", [(None, 3600), ("1", 1), ("14400", 14400)]
)
def test_timeout_default_and_boundaries_reach_prediction(
    monkeypatch, capsys, timeout, expected
):
    received = []

    def predict(args):
        received.append(args.timeout_seconds)
        return {"regions": [], "genes": []}

    monkeypatch.setattr(bounded, "predict", predict)
    args = PREDICT_ARGS + (["--timeout-seconds", timeout] if timeout else [])
    assert bounded.main(args) == 0
    assert received == [expected]
    assert capsys.readouterr().out == '{"regions": 0, "genes": 0}\n'


@pytest.mark.parametrize("timeout", ["0", "14401"])
def test_invalid_timeout_prevents_prediction(monkeypatch, capsys, timeout):
    def predict(args):
        pytest.fail("Invalid timeout must fail before prediction")

    monkeypatch.setattr(bounded, "predict", predict)
    with pytest.raises(SystemExit) as exc:
        bounded.main(PREDICT_ARGS + ["--timeout-seconds", timeout])
    assert exc.value.code == 2
    assert "--timeout-seconds: invalid choice" in capsys.readouterr().err

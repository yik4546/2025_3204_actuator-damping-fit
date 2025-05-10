% 파일 경로 목록
files = {'C:\Users\USER\actuator-fitting\data\1.csv', 'C:\Users\USER\actuator-fitting\data\2.csv', 'C:\Users\USER\actuator-fitting\data\3.csv', ...
        'C:\Users\USER\actuator-fitting\data\4.csv', 'C:\Users\USER\actuator-fitting\data\5.csv'};

% 초기 추정값과 최적화 옵션 설정
initial_c = 0.1;  % 감쇠 계수 c의 초기 추정값
options = optimset('Display', 'off', 'TolFun', 1e-6, 'TolX', 1e-6, ...
                   'MaxIter', 500, 'MaxFunEvals', 1000);

% 최적화된 파라미터와 매칭 퍼센지를 저장할 배열
params_all = NaN(length(files), 1);  % 실패한 경우 NaN으로 처리
match_all = zeros(length(files), 1);

% 각 파일에 대해 루프 실행
for i = 1:length(files)
    try
        % 실험 데이터 로드
        data = readmatrix(files{i});
        t_exp = data(:, 1);
        y_exp = data(:, 2);

        % 시간과 해당 데이터 정렬
        [t_exp, sortIdx] = sort(t_exp);
        y_exp = y_exp(sortIdx);
        y_exp = y_exp - y_exp(1);  % 실험 데이터의 초기값을 0으로 맞추기
        
        % 중복된 시간값이나 데이터를 제거
        t_exp = round(t_exp, 6);  % 시간 값을 소수점 6자리로 반올림
        [t_exp, unique_idx] = unique(t_exp, 'stable');  % 고유한 시간값만 선택
        y_exp = y_exp(unique_idx);  % 고유한 데이터만 남기기
        
        % 초기 높이와 길이 계산
        h0 = y_exp(1);           % 초기 높이
        L0 = h0 + 0.0228;        % 초기 높이에 기반한 길이

        % 모델 함수 정의
        model = @(c, t) simulate_actuator(c, t, h0, L0);

        % fminsearch를 사용하여 최적화 수행
        [c_opt, fval, exitflag] = fminsearch(@(c) sum((model(c, t_exp) - y_exp).^2), initial_c, options);
        
        % 최적화가 성공했는지 확인
        if exitflag <= 0
            fprintf('File %d: Optimization failed (exitflag = %d). Skipping...\n', i, exitflag);
            continue;  % 최적화 실패 시 해당 파일 건너뛰기
        end
        
        y_sim_opt = model(c_opt, t_exp);

        % 매칭 퍼센트 계산
        residuals = abs(y_exp - y_sim_opt);
        total_error = sum(residuals);
        total_data = sum(abs(y_exp));

        % 안전한 방식으로 매칭 퍼센트 계산
        if total_data > 0
            matching_percentage = (1 - total_error / total_data) * 100;
        else
            matching_percentage = 0;  % total_data가 0일 때는 0으로 설정
        end

        % 매칭 퍼센트가 임계값 이상일 경우 결과 저장
        if matching_percentage >= 0
            params_all(i) = c_opt;  % 감쇠 계수 (c) 저장
            match_all(i) = matching_percentage;

            % 피팅 결과를 플로팅
            figure;
            plot(t_exp, y_exp, 'ro-', 'DisplayName', '실험 데이터'); hold on;
            plot(t_exp, y_sim_opt, 'b-', 'DisplayName', '피팅된 시뮬레이션');
            title(sprintf('File %d Fitting Result - Matching: %.2f%%', i, matching_percentage));
            xlabel('시간 (s)');
            ylabel('높이 (m)');
            legend;
            grid on;

            fprintf('File %d: Matching Percentage = %.2f%%\n', i, matching_percentage);
        else
            fprintf('File %d: Matching below threshold (%.2f%%), skipped.\n', i, matching_percentage);
        end
    catch ME
        fprintf('File %d에서 오류 발생: %s\n', i, ME.message);
        if contains(ME.message, '샘플 점은 고유해야 합니다.')  % 오류 메시지가 해당하는 경우
            fprintf('File %d: 중복된 샘플 점으로 인해 피팅 실패. 건너뛰기...\n', i);
        else
            warning("File %d에서 피팅 실패: %s", i, ME.message);
        end
        continue;  % 오류 발생 시에도 파일 건너뛰기
    end
end

% 유효한 피팅된 파라미터들의 평균 계산
valid_params = params_all(~isnan(params_all));  % NaN 값 제거
if ~isempty(valid_params)
    avg_c = mean(valid_params);  % 유효한 감쇠 계수들의 평균 계산
    fprintf('\n유효한 피팅된 감쇠 계수의 평균:\n c = %.4f\n', avg_c);
else
    fprintf('\n유효한 데이터가 없으므로 피팅을 수행할 수 없습니다.\n');
end

% 액추에이터 동역학을 시뮬레이션하는 함수
function y = simulate_actuator(c, tvec, h, L)
    % 물리 상수
    mi = 0.00588;
    mo = 0.02285;
    ks = 455;
    alpha = 0.4937;
    g = 9.81;  % 중력 가속도 (m/s^2)

    % 초기 상태 (yi0 = h0로 설정)
    yi0 = h;  % 실험 데이터의 첫 번째 값 h0로 초기화
    yo0 = h + L;
    vi0 = 0;
    vo0 = 0;
    y0 = [yi0; yo0; vi0; vo0];

    % 시스템 행렬
    M = diag([mi, mo]);
    C = [c, -c; -c, c];
    K = [ks, -ks; -ks, ks];
    Fg = [-mi * g; -mo * g];

    % 시스템 방정식 (중력 포함)
    dydt = @(t, y) [y(3:4); M \ (-C * y(3:4) - K * y(1:2) + Fg)];

    % 이벤트 함수 정의 (yi = 0 도달 시 적분 중단)
    function [value, isterminal, direction] = groundEvent(t, y)
        value = y(1);           % yi = 0 (바닥과 충돌 감지)
        isterminal = 1;         % 적분 중단
        direction = -1;         % 위→아래로만 충돌 감지
    end

    options = odeset('Events', @groundEvent, 'RelTol', 1e-8, 'AbsTol', 1e-8);

    % 시뮬레이션 반복 수행
    t0 = tvec(1);  % 시뮬레이션 시간 시작
    tf = tvec(end);  % 시뮬레이션 시간 끝
    T = [];
    Y = [];

    while t0 < tf
        tspan = [t0, tf];

        [tSol, ySol, te, ye] = ode45(dydt, tspan, y0, options);

        T = [T; tSol];
        Y = [Y; ySol];

        if isempty(te)
            break;
        end

        % 충돌 시 반사 적용
        y0 = ye';
        y0(1) = 0;  % 충돌 후 y(1) = 0
        y0(3) = -alpha * y0(3);  % 반사 적용 (속도 방향 반전)
        t0 = te;
    end

    % 결과 보간 전에 중복된 시간 값을 처리합니다
    [T, unique_idx] = unique(T, 'stable');  % T에서 중복된 값 제거
    Y = Y(unique_idx, :);  % Y도 해당 인덱스를 통해 재정렬

    % tvec에 대해서도 중복값 처리
    [tvec, unique_tvec_idx] = unique(tvec, 'stable');  % tvec에서 중복된 값 제거

    % 결과 보간
    y = interp1(T, Y(:,1), tvec, 'linear', 'extrap');  % 보간 수행
end